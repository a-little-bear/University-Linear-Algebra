# Chapter 13: Matrix Functions

<div class="context-flow" markdown>

**Prerequisites**: Eigendecomposition (Ch10) · Jordan Form (Ch12) · Power Series

**Chapter Outline**: Definition of Matrix Power Series → Convergence → Matrix Exponential $e^A$ → Matrix Logarithm and Trigonometric Functions → Computation via Diagonalization → Computation via Jordan Form → Extension of Cauchy Integral Formula to Matrices → Applications (Solving Systems of Linear DEs)

**Extension**: Matrix exponential is the core operator in modern control theory and quantum mechanics, mapping algebraic addition to group multiplication.

</div>

Matrix functions extend the domain of a scalar function $f(z)$ to square matrices. This extension is not element-wise but maintains algebraic consistency. The most central tool is the **matrix exponential** $e^A$, which is the key to solving all linear continuous dynamical systems.

---

## 13.1 Definitions and Computation

!!! definition "Definition 13.1 (Matrix Power Series)"
    If a scalar function $f(z) = \sum_{k=0}^\infty a_k z^k$ converges in some disk, the matrix function is defined as:
    $$f(A) = \sum_{k=0}^\infty a_k A^k$$

!!! theorem "Theorem 13.3 (Diagonalization Method)"
    If $A = PDP^{-1}$, then $f(A) = P f(D) P^{-1} = P \operatorname{diag}(f(\lambda_i)) P^{-1}$.

---

## Exercises

1. **[Basic Calculation] If $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$, calculate $e^A$.**
   ??? success "Solution"
       $e^A = \begin{pmatrix} e^1 & 0 \\ 0 & e^2 \end{pmatrix}$. For diagonal matrices, simply apply the function to each diagonal element.

2. **[Exponential Properties] Prove: If $AB = BA$, then $e^{A+B} = e^A e^B$.**
   ??? success "Solution"
       Use the power series expansion and the binomial theorem. Since $A$ and $B$ commute, $(A+B)^k$ can be expanded like scalars. If they don't commute, this is generally false (requiring the BCH formula).

3. **[Nilpotent Matrix] Calculate $e^{At}$ where $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$.**
   ??? success "Solution"
       Since $A^2 = 0$, the series is $e^{At} = I + At + 0 = \begin{pmatrix} 1 & t \\ 0 & 1 \end{pmatrix}$.

4. **[Trace Identity] Prove $\det(e^A) = e^{\operatorname{tr}(A)}$.**
   ??? success "Solution"
       Let the eigenvalues of $A$ be $\lambda_i$. The eigenvalues of $e^A$ are $e^{\lambda_i}$.
       $\det(e^A) = \prod e^{\lambda_i} = e^{\sum \lambda_i} = e^{\operatorname{tr}(A)}$.

5. **[Trigonometric] If $A$ is a projection matrix ($A^2=A$), calculate $\sin(\pi A)$.**
   ??? success "Solution"
       $A^k = A$ for all $k \ge 1$.
       $\sin(\pi A) = \sum \frac{(-1)^k}{(2k+1)!} (\pi A)^{2k+1} = A \sum \frac{(-1)^k \pi^{2k+1}}{(2k+1)!} = A \sin(\pi) = 0$.

6. **[Jordan Block Function] Write the formula for $f(J_2(\lambda))$.**
   ??? success "Solution"
       $f \begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix} = \begin{pmatrix} f(\lambda) & f'(\lambda) \\ 0 & f(\lambda) \end{pmatrix}$. This reflects the relationship between the off-diagonal coupling and the derivative.

7. **[Invertibility] Is $e^A$ always invertible?**
   ??? success "Solution"
       Yes. Because its determinant $\det(e^A) = e^{\operatorname{tr}(A)}$ is never zero. Its inverse is $e^{-A}$.

8. **[Differentiation] Prove $\frac{d}{dt} e^{At} = A e^{At}$.**
   ??? success "Solution"
       Differentiating the series $\sum \frac{t^k A^k}{k!}$ term-by-term gives $\sum \frac{k t^{k-1} A^k}{k!} = A \sum \frac{t^{k-1} A^{k-1}}{(k-1)!} = A e^{At}$.

9. **[Rotation] Calculate $e^{Jt}$ where $J = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$.**
   ??? success "Solution"
       Since $J^2 = -I$, separating the series into even and odd terms gives the rotation matrix: $e^{Jt} = \begin{pmatrix} \cos t & -\sin t \\ \sin t & \cos t \end{pmatrix}$.

10. **[Application] How do you use matrix exponential to solve $\dot{x} = Ax, x(0)=x_0$?**
    ??? success "Solution"
        The solution is $x(t) = e^{At} x_0$. The matrix exponential maps the initial state to the state at any time $t$.

## Chapter Summary

Matrix functions are the advanced calculus of linear algebra:

1. **Structural Mapping**: They translate analytical properties of scalars perfectly into the operator space.
2. **Computational Core**: Diagonalization and Jordan chains are the standard paths for computing any matrix function.
3. **Dynamical Value**: The matrix exponential is the "time operator" for continuous systems, the ultimate tool for describing evolution.
