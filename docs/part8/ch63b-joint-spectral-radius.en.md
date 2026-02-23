# Chapter 63B: Joint Spectral Radius (JSR)

<div class="context-flow" markdown>

**Prerequisites**: Spectral Radius (Ch14) · Matrix Norms (Ch15) · Simultaneous Triangularization (Ch63A)

**Chapter Outline**: From Single Matrices to Matrix Sets → Definition of the Joint Spectral Radius (JSR) → Asymptotic Growth of Product Sequences → Key Identity: The Rota-Strang Identity → The Generalized Spectral Radius Formula → The Berger-Wang Theorem (Relation to standard spectral radius) → Extremality Principles and Computational Hardness → Applications: Stability of Switched Linear Systems, Smoothness of Wavelets, and Global Attractors in Discrete Dynamics

**Extension**: The JSR measures the "worst-case" evolution rate of a system; it extends spectral theory from a single trajectory to all possible combinations of paths, serving as the ultimate mathematical criterion for studying dynamical systems with uncertain switching logic.

</div>

In Ch14, we saw that the spectral radius $\rho(A)$ determines the convergence of $A^k$. However, in many modern applications (e.g., autonomous driving, wavelet construction), a system may choose, at each step, one matrix from a set $\{A_1, \ldots, A_m\}$ to execute. The **Joint Spectral Radius** (JSR) describes the growth rate under such a "worst-case" path. It reveals that even if every individual matrix is stable, their arbitrary combination can still lead to exponential explosion. This chapter explores this concept that is as computationally challenging as it is vital for stability.

---

## 63B.1 Definition of the Joint Spectral Radius

!!! definition "Definition 63B.1 (Joint Spectral Radius $\rho(\mathcal{F})$)"
    For a finite set of matrices $\mathcal{F} = \{A_1, \ldots, A_m\}$, the **Joint Spectral Radius** is defined as:
    $$\rho(\mathcal{F}) = \lim_{k \to \infty} \max_{A \in \mathcal{F}^k} \|A\|^{1/k}$$
    where $\mathcal{F}^k$ is the set of all possible products of length $k$ using elements from $\mathcal{F}$.
    **Intuition**: It is the maximum average scaling factor per step among all possible switching paths.

---

## 63B.2 Core Theorems

!!! theorem "Theorem 63B.1 (Rota-Strang Identity)"
    $$\rho(\mathcal{F}) = \inf_{\|\cdot\|} \max_{A \in \mathcal{F}} \|A\|$$
    where the infimum is taken over all operator norms. This means JSR can be viewed as the maximum single-step gain in some "optimal coordinate system."

!!! theorem "Theorem 63B.2 (Berger-Wang Theorem)"
    For a finite set $\mathcal{F}$, the joint spectral radius equals the supremum of the spectral radii of all finite products:
    $$\rho(\mathcal{F}) = \limsup_{k \to \infty} \max_{A \in \mathcal{F}^k} \rho(A)^{1/k}$$

---

## 63B.3 Stability Criterion

!!! technique "Stability of Switched Systems"
    A switched linear system $x_{k+1} = A_{\sigma_k} x_k$ is absolutely stable under arbitrary switching $\sigma_k$ iff:
    $$\rho(\{A_1, \ldots, A_m\}) < 1$$

---

## Exercises

**1. [Basics] Determine the JSR of the set $\mathcal{F} = \{ \begin{pmatrix} 0.5 & 0 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 0 & 0.5 \end{pmatrix} \}$.**

??? success "Solution"
    **Analysis:**
    1. Products of these matrices always contain zeros.
    2. Any product $A_{i_1} \cdots A_{i_k}$ has a norm that does not exceed $0.5^k$ (if all factors are the same) or 0 (if they are mixed).
    3. The eigenvalues are always 0.5 or 0.
    **Conclusion**: $\rho(\mathcal{F}) = 0.5$. Since the spectral radii are consistent and there is no constructive interference, switching does not increase growth.

**2. [Pitfall] Give an example where $\rho(A_1) < 1$ and $\rho(A_2) < 1$ but $\rho(\{A_1, A_2\}) > 1$.**

??? success "Solution"
    **Construction:**
    Let $A_1 = \begin{pmatrix} 0 & 2 \\ 0 & 0 \end{pmatrix}$ and $A_2 = \begin{pmatrix} 0 & 0 \\ 2 & 0 \end{pmatrix}$.
    1. Individually, $\rho(A_1) = 0$ and $\rho(A_2) = 0$. Any power is 0.
    2. Consider the product $A_1 A_2 = \begin{pmatrix} 4 & 0 \\ 0 & 0 \end{pmatrix}$.
    3. The eigenvalue of the product is 4.
    4. The modulus of the product sequence grows like $4^k$.
    **Conclusion**: $\rho(\{A_1, A_2\}) \ge \sqrt{4} = 2 > 1$. This proves that even if subsystems are stable, improper switching can lead to exponential instability.

**3. [Calculation] If $\mathcal{F} = \{A\}$ (a singleton), what does the JSR equal?**

??? success "Solution"
    **Conclusion: It equals the standard spectral radius $\rho(A)$.**
    This is exactly the definition of Gelfand's formula (Ch14). JSR is the natural generalization of Gelfand's formula to sets.

**4. [Property] Prove $\rho(\mathcal{F}) \le \max_{A \in \mathcal{F}} \|A\|$ for any consistent matrix norm.**

??? success "Solution"
    **Proof:**
    1. For any product $P = A_{i_1} \cdots A_{i_k}$, by sub-multiplicativity:
    2. $\|P\| \le \|A_{i_1}\| \cdots \|A_{i_k}\| \le (\max \|A\|)^k$.
    3. Taking the $k$-th root and the limit yields $\rho(\mathcal{F}) \le \max \|A\|$.

**5. [Hardness] What is the computational complexity of determining JSR?**

??? success "Solution"
    **Conclusion: Extremely difficult (NP-hard, even undecidable).**
    Determining whether $\rho(\mathcal{F}) \le 1$ has been proven undecidable in general. In practice, one uses Semidefinite Programming (SDP) to find quadratic Lyapunov functions to obtain upper bounds.

**6. [Wavelets] How does the JSR determine the smoothness of wavelet functions?**

??? success "Solution"
    A wavelet subdivision scheme is essentially a switched linear system. The regularity (differentiability) of the wavelet is directly determined by the JSR of the associated mask matrix set. A smaller JSR implies faster convergence of the subdivision and a smoother function.

**7. [Consistency] Prove: If the matrices in $\mathcal{F}$ are simultaneously triangularizable, then $\rho(\mathcal{F}) = \max \rho(A_i)$.**

??? success "Solution"
    **Reasoning:**
    If simultaneously triangularizable, the diagonal entries of any product are simply the products of the diagonal entries of the factors in the same basis. The largest possible eigenvalue must come from a combination of the largest eigenvalues of individual matrices. Since they commute (or are solvable), there is no cross-term reinforcement as seen in Problem 2.

**8. [Convergence] If $\rho(\mathcal{F}) = 0.99$, is the system necessarily stable?**

??? success "Solution"
    **Yes.**
    This means any product of length $k$ will eventually behave like $0.99^k$, converging to zero as $k$ goes to infinity.

**9. [Extremality] What is an "Extremal Cycle"?**

??? success "Solution"
    If there exists a finite-length product $A_{i_1} \cdots A_{i_k}$ whose spectral radius exactly matches the growth rate dictated by the JSR ($\rho(P)^{1/k} = \rho(\mathcal{F})$), it is called an extremal cycle. It represents the most dangerous switching mode of the system.

**10. [Application] Briefly state the significance of JSR in autonomous obstacle avoidance.**

??? success "Solution"
    Autonomous controllers switch laws depending on road conditions. To guarantee the vehicle remains controlled regardless of environment changes, engineers must ensure the JSR of the set of all possible control matrices is less than 1. This is the gold standard for robust stability in switched systems.

## Chapter Summary

The Joint Spectral Radius is the ultimate measure for uncertain dynamical systems:

1.  **Early Warning for Worst-case**: It reveals the boundaries of survival in dynamic environments, proving that local stability does not imply global safety without considering algebraic coherence between paths.
2.  **Unification of Norms and Spectra**: The Rota-Strang identity transforms complex asymptotic limits into a static search for an optimal norm, providing a theoretical pivot for stability assessment.
3.  **Scientific Limits of Computation**: The hardness of JSR marks the intersection of linear algebra and complexity theory, reminding us that when facing systems with infinite switching possibilities, we must rely on advanced optimization tools like convex relaxations.
