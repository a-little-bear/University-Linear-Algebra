# Chapter 14: Matrix Analysis

<div class="context-flow" markdown>

**Prerequisites**: Matrix Functions (Ch13) · Calculus Basics · Norms Basics (Ch15)

**Chapter Outline**: Matrix Sequences and Limits → Convergence of Operator Series → Matrix Calculus Introduction → Derivatives and Integrals → Relation between Spectral Radius and Convergence → Gelfand's Formula → Applications (Convergence of Iterative Methods)

**Extension**: Matrix analysis is the intersection of numerical analysis and continuous dynamics, bringing continuous tools from calculus into the discrete world of matrices.

</div>

Matrix analysis studies the analytical properties of matrices as variables. The core questions lie in defining "closeness" of matrices and handling the evolution of matrix sequences.

---

## 14.1 Sequences and Convergence

!!! definition "Definition 14.1 (Convergence of Matrix Sequences)"
    A matrix sequence $\{A_k\}$ converges to $A$ if for any matrix norm, $\lim_{k \to \infty} \|A_k - A\| = 0$.

!!! theorem "Theorem 14.3 (Gelfand's Formula)"
    For any matrix norm, the spectral radius $\rho(A)$ satisfies:
    $$\rho(A) = \lim_{k \to \infty} \|A^k\|^{1/k}$$

---

## Exercises

1. **[Convergence] Let $A = \begin{pmatrix} 0.5 & 1 \\ 0 & 0.5 \end{pmatrix}$. Does the sequence $A^k$ converge to the zero matrix? Explain.**
   ??? success "Solution"
       Yes. Both eigenvalues of $A$ are 0.5. Since the spectral radius $\rho(A) = 0.5 < 1$, the matrix power $A^k$ must converge to zero as $k \to \infty$.

2. **[Power Series] Given $A = \begin{pmatrix} 0.1 & 0 \\ 0 & 0.2 \end{pmatrix}$. Calculate the series $\sum_{k=0}^\infty A^k$.**
   ??? success "Solution"
       Since $\rho(A) = 0.2 < 1$, the series converges to $(I-A)^{-1}$.
       $I-A = \begin{pmatrix} 0.9 & 0 \\ 0 & 0.8 \end{pmatrix}$.
       Thus $\sum A^k = \begin{pmatrix} 1/0.9 & 0 \\ 0 & 1/0.8 \end{pmatrix} = \begin{pmatrix} 10/9 & 0 \\ 0 & 1.25 \end{pmatrix}$.

3. **[Spectral Radius] Provide an example of a matrix $A$ such that $\|A\|_2 > 1$ but $\rho(A) < 1$. What is the implication for dynamical systems?**
   ??? success "Solution"
       Take $A = \begin{pmatrix} 0 & 10 \\ 0 & 0 \end{pmatrix}$. The eigenvalues are 0, so $\rho(A)=0 < 1$. However, the norm $\|A\|_2 = 10 > 1$.
       This implies the system might experience large transient growth in the short term, but is ultimately stable in the long run.

4. **[Derivative] Calculate the derivative $A'(t)$ for $A(t) = \begin{pmatrix} t & t^2 \\ 1 & e^t \end{pmatrix}$.**
   ??? success "Solution"
       $A'(t) = \begin{pmatrix} 1 & 2t \\ 0 & e^t \end{pmatrix}$. The derivative of a matrix is the matrix of its component-wise derivatives.

5. **[Product Rule] Prove: $\frac{d}{dt}(A(t)B(t)) = A'(t)B(t) + A(t)B'(t)$.**
   ??? success "Solution"
       Expand $(A+\Delta A)(B+\Delta B) - AB$, drop high-order terms, and take the limit. Since matrix multiplication is non-commutative, the order of terms must be strictly maintained.

6. **[Integral] Calculate $\int_0^1 \begin{pmatrix} x & 1 \\ 0 & x^2 \end{pmatrix} dx$.**
   ??? success "Solution"
       $\begin{pmatrix} \int x & \int 1 \\ \int 0 & \int x^2 \end{pmatrix} = \begin{pmatrix} 1/2 & 1 \\ 0 & 1/3 \end{pmatrix}$.

7. **[Gelfand's App] If $\|A\| < 1$ for some norm, can you conclude $\rho(A) < 1$?**
   ??? success "Solution"
       Yes. Since the inequality $\rho(A) \le \|A\|$ holds for any induced norm, if the norm is less than 1, the spectral radius must be less than 1.

8. **[Inverse Derivative] Using $A(t)A^{-1}(t) = I$, prove $\frac{d}{dt}(A^{-1}(t)) = -A^{-1} A' A^{-1}$.**
   ??? success "Solution"
       Differentiate both sides: $A' A^{-1} + A (A^{-1})' = 0 \implies A (A^{-1})' = -A' A^{-1}$.
       Left-multiplying by $A^{-1}$ gives $(A^{-1})' = -A^{-1} A' A^{-1}$. Note $A'$ cannot be moved.

9. **[Trace Derivative] Prove $\frac{d}{dt} \operatorname{tr}(A(t)) = \operatorname{tr}(A'(t))$.**
   ??? success "Solution"
       Trace is a linear combination of entries; differentiation and summation operations commute.

10. **[Application] Why is the spectral radius of $B$ crucial in iterative methods $x_{k+1} = Bx_k + f$?**
    ??? success "Solution"
        The error evolution follows $e_{k+1} = B e_k \implies e_k = B^k e_0$. From matrix analysis, the error vanishes if and only if $\rho(B) < 1$. A smaller spectral radius implies faster convergence.

## Chapter Summary

Matrix analysis adds the dimension of time to static algebra:

1. **Convergence Criteria**: Spectral radius is the ultimate judge of invariant evolution.
2. **Analytical Tools**: Matrix calculus provides the symbolic engine for optimization and variational problems.
3. **Global Perspective**: Gelfand's formula perfectly unifies local norm estimates with global spectral attributes.
