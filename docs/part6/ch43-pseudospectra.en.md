# Chapter 43: Pseudospectra

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Matrix Norms (Ch15) · Resolvent (Ch13) · Matrix Stability (Ch36)

**Chapter Outline**: Limitation of Eigenvalues for Non-normal Matrices → Definition of $\epsilon$-Pseudospectrum → Resolvent Norm $\|(zI-A)^{-1}\|$ → Properties of Pseudospectra → Pseudospectral Boundary and Contours → Transient Growth vs. Asymptotic Decay → Kreiss Matrix Theorem → Numerical Computation

**Extension**: Pseudospectra explain why a system $\dot{x} = Ax$ can experience massive growth before eventually decaying, even if all eigenvalues are negative.

</div>

For **non-normal matrices** ($AA^* 
eq A^*A$), the spectrum $\sigma(A)$ can be extremely sensitive to perturbations. A small error in the entries can move the eigenvalues by a large amount. **Pseudospectra** provide a more robust description of an operator's behavior by identifying the regions of the complex plane where $(zI - A)$ is "nearly singular." This theory is essential for understanding transient phenomena in fluid dynamics, numerical analysis, and the behavior of iterative solvers.

---

## 43.1 Definitions and the Resolvent

!!! definition "Definition 43.1 ($\epsilon$-Pseudospectrum)"
    The $\epsilon$-pseudospectrum $\sigma_\epsilon(A)$ of a matrix $A$ is the set of $z \in \mathbb{C}$ satisfying any of the following equivalent conditions:
    1. $z$ is an eigenvalue of $A + E$ for some perturbation $\|E\| < \epsilon$.
    2. $\|(zI - A)^{-1}\| > 1/\epsilon$.
    3. $\sigma_{\min}(zI - A) < \epsilon$.

!!! theorem "Theorem 43.1 (Transient Growth)"
    Let $A$ be a Hurwitz stable matrix (all $\operatorname{Re}(\lambda) < 0$). If $\sigma_\epsilon(A)$ protrudes far into the right-half plane, the system $\dot{x} = Ax$ will exhibit **transient growth** before eventual decay:
    $$\sup_{t \ge 0} \|e^{At}\| \ge \frac{1}{\epsilon} \operatorname{dist}(\sigma_\epsilon(A), i\mathbb{R})$$

---

## Exercises

1. **[Fundamentals] Compute $\sigma_\epsilon(A)$ for a normal matrix $A$.**
   ??? success "Solution"
       For a normal matrix, the pseudospectrum is simply the union of open disks of radius $\epsilon$ centered at the eigenvalues: $\sigma_\epsilon(A) = \{z : \operatorname{dist}(z, \sigma(A)) < \epsilon\}$. Non-normality makes these sets much larger and non-spherical.

2. **[Resolvent] Relate the pseudospectrum to the resolvent operator $R(z, A) = (zI - A)^{-1}$.**
   ??? success "Solution"
       The pseudospectrum is the level set of the norm of the resolvent. Large resolvent norms indicate that $z$ is a "pseudo-eigenvalue," where the operator acts nearly singularly.

3. **[Transient Growth] Can a stable matrix have $\|e^{At}\| > 1$ for some $t > 0$?**
   ??? success "Solution"
       Yes, if the matrix is non-normal. Even if all eigenvalues satisfy $\operatorname{Re}(\lambda) < 0$, the non-orthogonal eigenvectors can interfere to produce massive transient growth before the asymptotic decay takes over. Pseudospectra quantify this risk.

4. **[Calculation] Let $A = \begin{pmatrix} -1 & 100 \ 0 & -1 \end{pmatrix}$. Where does $\sigma_{0.1}(A)$ lie?**
   ??? success "Solution"
       $A$ is highly non-normal. While $\sigma(A) = \{-1\}$, the pseudospectrum $\sigma_{0.1}(A)$ is a large region that extends into the right-half plane because $\|(zI-A)^{-1}\|$ is large for $z$ near $-1$ due to the large off-diagonal element.

5. **[Kreiss] State the Kreiss Matrix Theorem.**
   ??? success "Solution"
       $\sup_{t \ge 0} \|e^{At}\| \le e n \sup_{z \in \mathbb{C}^+} \operatorname{Re}(z) \|(zI-A)^{-1}\|$. This links the resolvent norm (pseudospectrum) to the maximum possible transient growth of the system.

6. **[Boundary] How does the boundary of $\sigma_\epsilon(A)$ relate to the singular values?**
   ??? success "Solution"
       The boundary $\partial \sigma_\epsilon(A)$ is the contour where the smallest singular value $\sigma_{\min}(zI - A)$ equals $\epsilon$.

7. **[Non-normality] Define the Henrici measure of non-normality.**
   ??? success "Solution"
       $v(A) = (\|A\|_F^2 - \sum |\lambda_i|^2)^{1/2}$. Larger $v(A)$ typically corresponds to larger and more distorted pseudospectra.

8. **[Numerical] How is the pseudospectrum computed in practice?**
   ??? success "Solution"
       By computing $\sigma_{\min}(zI - A)$ on a grid of points in the complex plane and plotting the contours. For large matrices, specialized algorithms like the Lanczos-based method for singular values are used.

9. **[Iterative Solvers] Why do pseudospectra predict the convergence of GMRES better than eigenvalues?**
   ??? success "Solution"
       For non-normal matrices, the convergence of GMRES is determined by the size of polynomials on the pseudospectrum, not just at the eigenvalues. If the pseudospectrum is large, GMRES will stagnate for many iterations before converging.

10. **[Spectral Mapping] Does $\sigma_\epsilon(f(A)) = f(\sigma_\epsilon(A))$?**
    ??? success "Solution"
        No. The Spectral Mapping Theorem does not hold for pseudospectra. The pseudospectrum of a function of $A$ can be much more complex than the image of the pseudospectrum under that function.

## Chapter Summary

This chapter establishes the robust theory of operator behavior under perturbation:

1. **Spectral Sensitivity**: Defined pseudospectra as the sets of possible eigenvalues under $\epsilon$-bounded perturbations.
2. **Resolvent Analysis**: Formulated the pseudospectrum as the level set of the resolvent norm, capturing near-singularity.
3. **Transient Dynamics**: Discovered the link between pseudospectral protrusion and the transient growth of the matrix exponential.
4. **Non-normal Metrics**: Positioned pseudospectra as the definitive tool for assessing the stability and convergence of non-normal operators.
