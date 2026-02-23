# Chapter 15: Norms and Perturbation Theory

<div class="context-flow" markdown>

**Prerequisites**: Inner Product Spaces (Ch8) · SVD (Ch11) · Matrix Analysis (Ch14)

**Chapter Outline**: Vector Norms ($L_1, L_2, L_\infty$) → Induced Matrix Norms → Operator and Spectral Norms → Frobenius Norm → Equivalence of Norms → Condition Number $\kappa(A)$ → Perturbation Theory (Linear Systems and Eigenvalues) → Bauer-Fike Theorem

**Extension**: Norms are rulers for measuring the "size" of mathematical objects, while perturbation theory studies how results fluctuate when real-world noise interferes with the input.

</div>

In pure mathematics, we talk about exact solutions, but in numerical linear algebra, we talk about error. **Norms** quantify the size of the error, and the **condition number** reveals the system's sensitivity to that error.

---

## 15.1 Core Definitions and Inequalities

!!! definition "Definition 15.1 (Vector $p$-norm)"
    For a vector $x$, its $p$-norm is defined as:
    $$\|x\|_p = \left( \sum |x_i|^p \right)^{1/p}$$
    Common cases are $p=1, 2, \infty$.

!!! theorem "Theorem 15.3 (Bauer-Fike Theorem)"
    If $A$ is diagonalizable ($A = VDV^{-1}$) and $\mu$ is an eigenvalue of $A+E$, then there exists an eigenvalue $\lambda$ of $A$ such that:
    $$|\mu - \lambda| \le \kappa_p(V) \|E\|_p$$

---

## Exercises

1. **[Vector Norms] Calculate the $L_1, L_2, L_\infty$ norms of $x = (3, -4)^T$.**
   ??? success "Solution"
       - $\|x\|_1 = |3| + |-4| = 7$.
       - $\|x\|_2 = \sqrt{3^2 + (-4)^2} = 5$.
       - $\|x\|_\infty = \max(|3|, |-4|) = 4$.

2. **[Frobenius] Calculate the Frobenius norm of $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$.**
   ??? success "Solution"
       $\|A\|_F = \sqrt{1^2 + 2^2 + 3^2 + 4^2} = \sqrt{30} \approx 4.47$.

3. **[Spectral Norm] What is the 2-norm (spectral norm) of $A = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$?**
   ??? success "Solution"
       For a diagonal matrix, the 2-norm is the maximum absolute value of the diagonal entries. Thus $\|A\|_2 = 3$.

4. **[Equivalence] Prove: In finite-dimensional spaces, any two norms $\|\cdot\|_a$ and $\|\cdot\|_b$ are equivalent.**
   ??? success "Solution"
       Since the unit sphere is compact in one norm, the other norm (as a continuous function) must have a maximum and minimum on this set. This guarantees constants $C_1, C_2$ such that $C_1 \|x\|_a \le \|x\|_b \le C_2 \|x\|_a$.

5. **[Condition Number] If $A = \begin{pmatrix} 1 & 0 \\ 0 & 0.01 \end{pmatrix}$, calculate its condition number with respect to the 2-norm.**
   ??? success "Solution"
       $\kappa_2(A) = \|A\|_2 \|A^{-1}\|_2 = 1 \cdot (1/0.01) = 100$. This indicates that the matrix amplifies errors 100 times during inversion.

6. **[Submultiplicativity] Prove that induced matrix norms satisfy $\|AB\| \le \|A\| \|B\|$.**
   ??? success "Solution"
       $\|ABx\| \le \|A\| \|Bx\| \le \|A\| (\|B\| \|x\|)$.
       By definition, $\|AB\| = \max \frac{\|ABx\|}{\|x\|} \le \|A\| \|B\|$.

7. **[Perturbation Bound] If input $b$ in $Ax=b$ is perturbed by $\Delta b$, what is the relative error upper bound for the solution?**
   ??? success "Solution"
       $\frac{\|\Delta x\|}{\|x\|} \le \kappa(A) \frac{\|\Delta b\|}{\|b\|}$. The condition number is the amplification factor for error propagation.

8. **[Eigenvalue Sensitivity] Why are eigenvalues of normal matrices (like symmetric matrices) more robust than non-normal ones?**
   ??? success "Solution"
       For normal matrices, the diagonalizing matrix $V$ can be chosen as unitary, making $\kappa_2(V) = 1$. The Bauer-Fike Theorem simplifies to $|\mu-\lambda| \le \|E\|_2$, meaning the eigenvalue shift is bounded by the perturbation size.

9. **[Calculation] Find the 1-norm (max absolute column sum) of $\begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       Column 1 sum: $1+0=1$. Column 2 sum: $2+1=3$.
       Thus $\|A\|_1 = \max(1, 3) = 3$.

10. **[Application] Why do we prefer unitary transformations (like Householder) in numerical methods?**
    ??? success "Solution"
        Because unitary transformations have a spectral norm of 1 and a condition number of 1. They do not amplify rounding errors, ensuring numerical stability of the algorithm.

## Chapter Summary

Norms and perturbation theory are the red lines of computational mathematics:

1. **Size Measurement**: Norms translate matrix properties into comparable numerical values.
2. **Stability Determination**: Condition number is the only barometer for algorithmic reliability.
3. **Error Control**: Perturbation theory establishes the valid boundaries of linear algebra calculations in a precision-limited real world.
