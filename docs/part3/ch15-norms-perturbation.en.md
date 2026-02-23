# Chapter 15: Matrix Norms and Perturbation Theory

<div class="context-flow" markdown>

**Prerequisites**: Matrix Analysis (Ch14) · Eigenvalues (Ch6) · Singular Values (Ch11) · Vector Spaces (Ch4)

**Chapter Outline**: Vector Norms ($L_1, L_2, L_\infty$) → Induced Matrix Norms → Submultiplicativity → Spectral Norm vs. Frobenius Norm → Schatten Norms → Condition Number $\kappa(A)$ → Perturbation of Linear Systems → Error Bounds → Backward Stability Overview

**Extension**: Matrix norms quantify the "size" of an operator; the condition number is the definitive measure of how much rounding errors or noise in the data will be magnified in the solution.

</div>

Matrix norms extend the concept of length from vectors to operators. An induced matrix norm measures the maximum possible magnification a matrix can apply to a vector. This magnification factor is the core of **perturbation theory**: it allows us to estimate how errors in the input $b$ or the matrix $A$ affect the solution $x$ in $Ax=b$. The **condition number** $\kappa(A)$ is the critical sensitivity index that determines the numerical feasibility of solving a linear system.

---

## 15.1 Definitions and Structural Properties

!!! definition "Definition 15.1 (Induced Matrix Norm)"
    The matrix norm induced by a vector norm $\|\cdot\|$ is:
    $$\|A\| = \sup_{x \neq 0} \frac{\|Ax\|}{\|x\|}$$
    Common induced norms include the $L_1$ (max column sum), $L_\infty$ (max row sum), and $L_2$ (spectral norm, $\|A\|_2 = \sigma_{\max}(A)$).

!!! theorem "Theorem 15.1 (Sensitivity Bound)"
    In the system $Ax = b$, if the input $b$ is perturbed by $\delta b$, the relative change in the solution satisfies:
    $$\frac{\|\delta x\|}{\|x\|} \le \kappa(A) \frac{\|\delta b\|}{\|b\|}, \quad \text{where } \kappa(A) = \|A\| \|A^{-1}\|$$

---

## Exercises

1. **[Fundamentals] Compute $\|A\|_1$ and $\|A\|_\infty$ for $A = \begin{pmatrix} 1 & -5 \\ 2 & 3 \end{pmatrix}$.**
   ??? success "Solution"
       $\|A\|_1 = \max(1+2, |-5|+3) = 8$ (column sums). $\|A\|_\infty = \max(1+|-5|, 2+3) = 6$ (row sums).

2. **[Spectral Norm] Find $\|A\|_2$ for $A = \begin{pmatrix} 3 & 0 \\ 0 & -2 \end{pmatrix}$.**
   ??? success "Solution"
       $\|A\|_2$ is the largest singular value. $\sigma = \{3, 2\}$. Thus $\|A\|_2 = 3$.

3. **[Frobenius] Define the Frobenius norm and relate it to the trace.**
   ??? success "Solution"
       $\|A\|_F = \sqrt{\sum a_{ij}^2} = \sqrt{\operatorname{tr}(A^T A)}$. It is the Euclidean norm of the matrix viewed as a vector in $\mathbb{R}^{n^2}$.

4. **[Condition Number] Calculate $\kappa_2(A)$ for $A = \begin{pmatrix} 10 & 0 \\ 0 & 0.1 \end{pmatrix}$.**
   ??? success "Solution"
       $\|A\|_2 = 10$, $\|A^{-1}\|_2 = 1/0.1 = 10$. $\kappa_2(A) = 10 \times 10 = 100$.

5. **[Submultiplicativity] Prove $\|AB\| \le \|A\| \|B\|$ for induced norms.**
   ??? success "Solution"
       $\|ABx\| \le \|A\| \|Bx\| \le \|A\| (\|B\| \|x\|) = (\|A\| \|B\|) \|x\|$. Taking the supremum over unit $x$ gives the result. This property is vital for analyzing error propagation.

6. **[Unitary Invariance] Why is the Frobenius norm called unitarily invariant?**
   ??? success "Solution"
       Because $\|UAV\|_F = \|A\|_F$ for any unitary $U, V$. This means the "size" of the matrix is independent of the orthonormal basis used.

7. **[Consistency] Show that $\|Av\| \le \|A\| \|v\|$ for any consistent norm pair.**
   ??? success "Solution"
       This is the defining property of induced norms. It allows us to treat matrices as bounded operators on vector spaces.

8. **[Singular Matrices] What happens to $\kappa(A)$ as $A$ approaches a singular matrix?**
   ??? success "Solution"
       As $A$ becomes singular, $\sigma_{\min} \to 0$, so $\|A^{-1}\| \to \infty$. The condition number $\kappa(A)$ explodes to infinity, reflecting extreme sensitivity.

9. **[Eigenvalues] Is $\|A\|_2 = \rho(A)$ for all matrices?**
   ??? success "Solution"
       No. Only for normal matrices ($AA^* = A^*A$). For non-normal matrices, $\|A\|_2$ can be much larger than $\rho(A)$.

10. **[Perturbation] If $\kappa(A) = 10^6$ and your data has 8 digits of precision, how many digits can you trust in the solution $x$?**
    ??? success "Solution"
        Approximately $8 - \log_{10}(10^6) = 2$ digits. The condition number acts as a "precision destroyer" in linear solvers.

## Chapter Summary

This chapter explores the metric sensitivity of linear operators:

1. **Magnification Metrics**: Defined induced and entry-wise norms to quantify the scale of matrix actions.
2. **Spectral Connection**: Established the $L_2$ norm as the maximum singular value, linking geometry to singular spectra.
3. **Sensitivity Calculus**: Formulated the condition number as the definitive multiplier for error propagation.
4. **Numerical Feasibility**: Linked the geometric properties of the operator to the precision limits of computational algorithms.
