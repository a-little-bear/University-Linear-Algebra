# Chapter 16: Positive Definite Matrices

<div class="context-flow" markdown>

**Prerequisites**: Inner Product Spaces (Ch8) · Quadratic Forms (Ch9) · Eigenvalues (Ch6)

**Chapter Outline**: Definition of Positive Definite and Semi-definite Matrices → Criteria (Eigenvalues, Principal Minors, Cholesky) → Properties (Determinant, Inverse, Trace) → Matrix Square Root $\sqrt{A}$ → Generalized Eigenvalue Problems → Extremal Properties and Rayleigh Quotient → Applications (Energy Functions, Vibration Analysis)

**Extension**: Positive definite matrices are the "positive numbers" of linear spaces; they ensure the convexity of energy functions and are algebraic prerequisites for the stability of physical systems.

</div>

Positive definite matrices form a well-behaved cone in the space of matrices. They provide positive "feedback" in all directions, ensuring system stability and the existence of optimal solutions.

---

## 16.1 Definitions and Core Criteria

!!! definition "Definition 16.1 (Positive Definite Matrix)"
    A symmetric matrix $A$ is positive definite (denoted $A \succ 0$) if for all non-zero vectors $x$:
    $$x^T A x > 0$$

!!! theorem "Theorem 16.3 (Equivalent Characterizations)"
    For a symmetric matrix $A$, the following are equivalent:
    1. $A$ is positive definite.
    2. All eigenvalues of $A$ are strictly greater than 0.
    3. All leading principal minors of $A$ are positive.
    4. There exists an invertible matrix $L$ such that $A = LL^T$.

---

## Exercises

1. **[Criteria] Is $A = \begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix}$ positive definite?**
   ??? success "Solution"
       Yes. Leading principal minors: $D_1 = 1 > 0$, $D_2 = \det A = 5-4 = 1 > 0$. Since all are positive, $A$ is positive definite.

2. **[Eigenvalues] If the eigenvalues of $A$ are $0, 1, 2$, is it positive definite or semi-definite?**
   ??? success "Solution"
       It is **positive semi-definite** ($A \succeq 0$). Positive definiteness requires eigenvalues to be strictly greater than 0; because of the 0 eigenvalue, it is not positive definite.

3. **[Property] If $A, B \succ 0$, prove $A+B \succ 0$.**
   ??? success "Solution"
       For any non-zero $x$, $x^T(A+B)x = x^T Ax + x^T Bx$. Since $x^T Ax > 0$ and $x^T Bx > 0$, their sum must be positive.

4. **[Inverse] If $A \succ 0$, prove $A^{-1}$ exists and $A^{-1} \succ 0$.**
   ??? success "Solution"
       Since eigenvalues $\lambda_i > 0$, $\det A = \prod \lambda_i > 0$, so $A$ is invertible. The eigenvalues of $A^{-1}$ are $1/\lambda_i$. Since $1/\lambda_i > 0$, $A^{-1}$ is also positive definite.

5. **[Square Root] Calculate the positive definite square root $A^{1/2}$ for $A = \begin{pmatrix} 4 & 0 \\ 0 & 9 \end{pmatrix}$.**
   ??? success "Solution"
       $A^{1/2} = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$.

6. **[Gram Matrix] Prove: For any matrix $M$ with full column rank, $A = M^T M$ is positive definite.**
   ??? success "Solution"
       $x^T A x = x^T M^T M x = \|Mx\|^2$. Since $M$ has full column rank and $x \neq 0$, then $Mx \neq 0$, so $\|Mx\|^2 > 0$.

7. **[Rayleigh Quotient] Find the range of values for $\frac{x^T Ax}{x^T x}$ where $A = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       The range is $[\lambda_{\min}, \lambda_{\max}] = [1, 2]$.

8. **[Determinant] Prove: The determinant of a positive definite matrix is less than or equal to the product of its diagonal entries.**
   ??? success "Solution"
       This is a consequence of Hadamard's Inequality: $\det A \le \prod a_{ii}$, with equality iff $A$ is diagonal.

9. **[Schur Complement] Determine the condition for $a$ such that $\begin{pmatrix} 1 & 1 \\ 1 & a \end{pmatrix}$ is PD.**
   ??? success "Solution"
       Principal minors: $1 > 0$ and $a-1 > 0 \implies a > 1$.

10. **[Application] Why does a positive definite Hessian matrix guarantee convergence to a local minimum in gradient descent?**
    ??? success "Solution"
        A positive definite Hessian implies the objective function is locally convex (bowl-shaped) at the stationary point, ensuring the point is a unique local minimum rather than a maximum or saddle point.

## Chapter Summary

Positive definite matrices are the most "benign" objects in linear algebra:

1. **Spectral Positivity**: All eigenvalues being positive defines the expansionary property of the matrix.
2. **Energy Geometry**: As an energy functional, the convexity of $x^T Ax$ is fully determined by positive definiteness.
3. **Decomposition Advantage**: Cholesky decomposition is the most powerful tool for simplifying numerical calculations using positive definiteness.
