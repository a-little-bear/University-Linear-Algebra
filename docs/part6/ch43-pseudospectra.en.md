# Chapter 43: Pseudospectra and Non-normal Matrix Analysis

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Matrix Norms and Perturbations (Ch15) · Matrix Functions (Ch13)

**Chapter Outline**: The Crisis of Non-normality (Deceptive Eigenvalues) → Definition of $\epsilon$-Pseudospectra → Equivalent Characterizations (Resolvent Norm, Singular Values, Perturbed Spectrum) → Measuring Departure from Normality (Henrici’s Number) → Geometry of Pseudospectra in the Complex Plane → Transient Growth in Non-normal Systems → The Kreiss Matrix Theorem → Applications: Hydrodynamic Stability, Convergence of Numerical Algorithms, and Laser Resonator Analysis

**Extension**: Pseudospectra theory is a "generalized spectral theory" for highly non-normal operators; it reveals that when a matrix deviates significantly from being normal, traditional eigenvalue analysis fails. The dynamical behavior of the system is instead governed by the pseudospectrum—a larger set—providing the key to understanding "pseudo-stability" in physical systems.

</div>

For normal matrices, eigenvalues perfectly dictate the stability of a system. However, for highly **non-normal matrices** ($AA^* \neq A^*A$), eigenvalues can be extremely misleading. Even if all eigenvalues lie in the left half-plane, the system may still undergo massive transient growth. **Pseudospectra** provide a deeper description of these "ill-conditioned" operators by studying the norm of the resolvent. This chapter uncovers the hidden regions of potential instability lurking behind the eigenvalues.

---

## 43.1 Definition of $\epsilon$-Pseudospectra

!!! definition "Definition 43.1 ($\epsilon$-Pseudospectrum)"
    For $\epsilon > 0$, the **$\epsilon$-pseudospectrum** $\sigma_\epsilon(A)$ of a matrix $A$ is the set in the complex plane satisfying the following equivalent conditions:
    1.  **Perturbation**: There exists a matrix $E$ with $\|E\| < \epsilon$ such that $z \in \sigma(A+E)$.
    2.  **Resolvent**: The norm of its resolvent satisfies $\|(zI - A)^{-1}\| > 1/\epsilon$.
    3.  **Singular Value**: $s_{\min}(zI - A) < \epsilon$.

---

## 43.2 Non-normality and Transient Growth

!!! definition "Definition 43.2 (Departure from Normality)"
    The non-normality of matrix $A$ can be measured by **Henrici’s Number**: $\Delta(A) = \sqrt{\|A\|_F^2 - \|\Lambda\|_F^2}$.
    The stronger the non-normality, the further the pseudospectral region expands relative to the convex hull of the eigenvalues.

!!! theorem "Theorem 43.1 (Kreiss Matrix Theorem)"
    The peak of transient growth in the system $\|e^{At}\|$ is proportional to how far the pseudospectrum protrudes into the right half-plane. This explains why "theoretically stable" fluid flows can suddenly transition to turbulence.

---

## 43.3 Geometry of Pseudospectra

!!! technique "Visualization: Contour Plots"
    Pseudospectra are typically visualized by plotting the contours of $f(z) = \log_{10} s_{\min}(zI - A)$.
    - For **normal matrices**: Contours are perfect circles centered at the eigenvalues.
    - For **non-normal matrices**: Contours are heavily distorted and merge into large "danger zones."

---

## Exercises

**1. [Basics] Prove that for a normal matrix, the $\epsilon$-pseudospectrum is simply the $\epsilon$-neighborhood of its eigenvalues.**

??? success "Solution"
    **Proof:**
    1. For a normal matrix, $\|(zI - A)^{-1}\|_2 = 1 / \operatorname{dist}(z, \sigma(A))$.
    2. By definition, $z \in \sigma_\epsilon(A) \iff \|(zI - A)^{-1}\| > 1/\epsilon$.
    3. Substituting yields: $1 / \operatorname{dist}(z, \sigma(A)) > 1/\epsilon \iff \operatorname{dist}(z, \sigma(A)) < \epsilon$.
    **Conclusion**: The pseudospectrum of a normal matrix does not "spill over"; stability is entirely determined by eigenvalues.

**2. [Calculation] Estimate the $\epsilon=1$ pseudospectrum for $A = \begin{pmatrix} 0 & 100 \\ 0 & 0 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. Calculate $s_{\min}(zI - A) = s_{\min} \begin{pmatrix} z & -100 \\ 0 & z \end{pmatrix}$.
    2. The eigenvalues are $0, 0$.
    3. When $z$ is small, the matrix is near rank-1. The smallest singular value is approximately $|z|^2/100$.
    4. Setting $|z|^2/100 < 1 \implies |z| < 10$.
    **Conclusion**: Despite having zero as the only eigenvalue, the $\epsilon=1$ pseudospectrum is a disk with radius $\approx 10$. Tiny perturbations can shift the eigenvalues by as much as 10.

**3. [Property] Prove that $\sigma_\epsilon(A)$ always contains the spectrum $\sigma(A)$.**

??? success "Solution"
    **Proof:**
    1. When $z \in \sigma(A)$, the matrix $zI - A$ is singular.
    2. Its resolvent norm $\|(zI - A)^{-1}\|$ becomes infinite.
    3. For any $\epsilon > 0$, $\infty > 1/\epsilon$.
    **Conclusion**: Eigenvalues are the "source points" of the pseudospectral sets.

**4. [Transient Growth] Why doesn't having all eigenvalues in the left half-plane guarantee that $\|e^{At}\|$ won't increase before decaying?**

??? success "Solution"
    **Algebraic Reason:**
    1. The initial derivative of $\|e^{At}\|$ (bounded by the numerical range) can be positive.
    2. If $A$ is non-normal, the eigenvector basis is non-orthogonal.
    3. Components can interfere destructively initially. When this balance is disrupted, the system's energy can temporarily spike.
    4. Only in the limit $t \to \infty$ does the exponential decay of eigenvalues dominate.

**5. [Numerical] Briefly describe how to efficiently compute pseudospectra for large matrices.**

??? success "Solution"
    Since computing singular values across the entire plane is expensive:
    1. **Lanczos Projection**: Project the large matrix onto a smaller Krylov subspace.
    2. **Path-following Methods**: Use predictor-corrector methods to trace the contours of the resolvent norm.

**6. [Inclusion] If $\epsilon_1 < \epsilon_2$, prove $\sigma_{\epsilon_1}(A) \subset \sigma_{\epsilon_2}(A)$.**

??? success "Solution"
    **Proof:**
    If $z \in \sigma_{\epsilon_1}(A)$, then $\|(zI - A)^{-1}\| > 1/\epsilon_1$.
    Since $\epsilon_1 < \epsilon_2$, we have $1/\epsilon_1 > 1/\epsilon_2$.
    Thus $\|(zI - A)^{-1}\| > 1/\epsilon_2$, which satisfies $z \in \sigma_{\epsilon_2}(A)$.

**7. [Application] What role do pseudospectra play in numerical convergence analysis?**

??? success "Solution"
    **Significance:**
    For certain iterative operators (like non-symmetric preconditioners), even if the spectral radius $\rho(B) < 1$, strong pseudospectral effects due to non-normality can cause iterative residuals to explode initially, making the algorithm fail in practice. Pseudospectra analysis predicts this "early divergence" phenomenon.

**8. [Unitary Equivalence] Prove $\sigma_\epsilon(U^* A U) = \sigma_\epsilon(A)$ for unitary $U$.**

??? success "Solution"
    **Proof:**
    1. $\|(zI - U^*AU)^{-1}\| = \|(U^*(zI - A)U)^{-1}\| = \|U^*(zI - A)^{-1}U\|$.
    2. Since unitary transformations preserve the operator norm, this equals $\|(zI - A)^{-1}\|$.
    **Conclusion**: Pseudospectra are invariant under rotation of the coordinate system.

**9. [Henrici] Calculate Henrici’s non-normality measure $\Delta(A)$ for $A = \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. $\|A\|_F^2 = 1^2 + 2^2 + 0^2 + 3^2 = 14$.
    2. Eigenvalues are $1, 3$. $\|\Lambda\|_F^2 = 1^2 + 3^2 = 10$.
    3. $\Delta(A) = \sqrt{14 - 10} = 2$.
    **Conclusion**: The non-normality arises entirely from the off-diagonal entry 2.

**10. [Physics] Why was pseudospectra theory a breakthrough in fluid turbulence studies?**

??? success "Solution"
    **Explanation:**
    Classical stability theory predicted laminar flow should be stable at certain Reynolds numbers based on eigenvalues. Yet experiments showed transitions to turbulence much earlier. Pseudospectra theory proved that linearized fluid operators are **highly non-normal**. Huge expansions of the pseudospectrum into the right half-plane indicated that tiny disturbances would trigger massive transient amplification, the true algebraic culprit behind non-linear instability and turbulence.

## Chapter Summary

Pseudospectra theory is the modern diagnostic for the essence of instability:

1.  **Limits of Eigenvalues**: It proves that for highly coupled, non-symmetric systems, discrete spectral points cannot represent the full behavior of an operator; the energy response of the "near-spectrum" must be considered.
2.  **Power of the Resolvent**: By contouring the resolvent norm, pseudospectra transform abstract perturbation theory into a geometric map on the complex plane, establishing a new yardstick for numerical robustness.
3.  **Insight into Transients**: This chapter extends static canonical form analysis to dynamic energy analysis, revealing the "early-stage disasters" that can occur in non-normal systems before they asymptotically settle into stability.
