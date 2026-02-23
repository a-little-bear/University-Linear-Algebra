# Chapter 13: Matrix Functions

<div class="context-flow" markdown>

**Prerequisites**: Matrix Power (Ch2) · Jordan Form (Ch12) · Taylor Series · Diagonalization (Ch6)

**Chapter Outline**: Definition of $f(A)$ via Power Series → Definition via Jordan Form → Matrix Exponential $e^{At}$ → Matrix Logarithm $\ln(A)$ → Square Root of a Matrix $\sqrt{A}$ → Properties of Matrix Functions → Solving $\dot{x} = Ax$ → Commutativity and Matrix Functions → Cauchy Integral Formula for Matrices

**Extension**: Matrix functions generalize scalar calculus to operators; the matrix exponential is the fundamental solution to systems of linear differential equations.

</div>

Matrix functions allow us to apply standard mathematical functions—like $e^x, \sin x, \ln x$—directly to matrices. Instead of operating on individual entries, we operate on the operator as a whole. For a diagonalizable matrix $A = PDP^{-1}$, the definition is straightforward: $f(A) = P f(D) P^{-1}$. For general matrices, we use power series or the Jordan Canonical Form. The **matrix exponential** $e^{At}$ is the most critical function in this chapter, as it provides the explicit solution to the state-space equations of dynamical systems.

---

## 13.1 Defining Functions of Matrices

!!! definition "Definition 13.1 (Matrix Function via Series)"
    If $f(z)$ has a Taylor series $\sum a_k z^k$, the matrix function $f(A)$ is defined as:
    $$f(A) = \sum_{k=0}^\infty a_k A^k$$
    provided the series converges for all eigenvalues of $A$.

!!! theorem "Theorem 13.1 (Spectral Mapping Theorem)"
    The eigenvalues of $f(A)$ are exactly $\{f(\lambda_1), \dots, f(\lambda_n)\}$, where $\{\lambda_i\}$ are the eigenvalues of $A$.

---

## Exercises

1. **[Fundamentals] Compute $e^A$ for $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$.**
   ??? success "Solution"
       For a diagonal matrix, $e^A = \begin{pmatrix} e^1 & 0 \\ 0 & e^2 \end{pmatrix}$.

2. **[Series] Use the series definition to find $e^A$ for $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$.**
   ??? success "Solution"
       $A^2 = 0$. $e^A = I + A + \frac{1}{2!}A^2 + \dots = I + A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$.

3. **[Diagonalization] Compute $\sin A$ for $A = P \operatorname{diag}(\pi, \pi/2) P^{-1}$.**
   ??? success "Solution"
       $\sin A = P \operatorname{diag}(\sin \pi, \sin \pi/2) P^{-1} = P \operatorname{diag}(0, 1) P^{-1}$.

4. **[Logarithm] When does the matrix logarithm $\ln A$ exist?**
   ??? success "Solution"
       A real matrix has a real logarithm iff it is non-singular and every Jordan block associated with a negative eigenvalue occurs an even number of times (ensuring complex conjugate pairings). Generally, $A$ must be non-singular.

5. **[Differential Equations] Solve $\dot{x} = Ax$ with $x(0) = x_0$.**
   ??? success "Solution"
       The unique solution is $x(t) = e^{At} x_0$. The matrix exponential acts as the "propagator" of the system state.

6. **[Commutativity] Is $e^{A+B} = e^A e^B$ always true?**
   ??? success "Solution"
       No. It is true if and only if $A$ and $B$ commute ($AB = BA$).

7. **[Square Root] Find a square root of $A = \begin{pmatrix} 4 & 0 \\ 0 & 9 \end{pmatrix}$.**
   ??? success "Solution"
       One solution is $\operatorname{diag}(2, 3)$. Note that matrix square roots are not unique; $\operatorname{diag}(-2, 3)$ is also a square root.

8. **[Determinant] Prove the identity $\det(e^A) = e^{\operatorname{tr}(A)}$.**
   ??? success "Solution"
       $\det(e^A) = \prod e^{\lambda_i} = e^{\sum \lambda_i} = e^{\operatorname{tr}(A)}$. This link between the trace and determinant is fundamental in Lie theory.

9. **[Jordan Form] Describe how to compute $f(J_k(\lambda))$ for a Jordan block.**
   ??? success "Solution"
       $f(J_k(\lambda))$ is an upper triangular matrix where the diagonal is $f(\lambda)$, the first super-diagonal is $f'(\lambda)$, the second is $f''(\lambda)/2!$, and so on.

10. **[Inversion] Express $A^{-1}$ as a matrix function.**
    ??? success "Solution"
        $f(z) = 1/z$. For non-singular $A$, $A^{-1}$ can be computed via the Jordan form or power series (if within the radius of convergence).

## Chapter Summary

This chapter extends the domain of calculus to linear operators:

1. **Analytical Mapping**: Defined matrix functions through Taylor series and spectral decompositions.
2. **Spectral Propagation**: Established the Spectral Mapping Theorem as the rule for eigenvalue transformation.
3. **Dynamic Solution**: Developed the matrix exponential as the definitive analytical tool for linear differential systems.
4. **Algebraic Identity**: Leveraged the trace-determinant relation to link matrix functions to the core invariants of linear algebra.
