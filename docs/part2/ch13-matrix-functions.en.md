# Chapter 13: Matrix Functions

<div class="context-flow" markdown>

**Prerequisites**: Jordan Canonical Form (Ch12) · Matrix Analysis Basics (Ch14) · Calculus Series Theory

**Chapter Outline**: From Scalar to Matrix Functions → Power Series Definition → Jordan Form Definition Method → Sylvester-Lagrange Interpolation Method → Core Functions: Matrix Exponential ($e^A$), Matrix Logarithm ($\log A$), and Trigonometric Functions → Computation Techniques → Identities and Properties → Applications in Systems of ODEs

**Extension**: Matrix functions elevate arithmetic algebra to analytical algebra; they are the central mathematical language for control theory, evolution operators in quantum mechanics, and continuous dynamical systems.

</div>

Matrix functions study how to apply scalar functions (such as $e^x, \sin x, \log x$) to matrix variables. This is not a simple element-wise operation but an operator computation that preserves the algebraic structure of the matrix. Matrix functions bridge the gap between discrete matrix algebra and the continuous physical world.

---

## 13.1 Methods of Definition

!!! definition "Definition 13.1 (Definition via Jordan Form)"
    Let $A = P J P^{-1}$, where $J = \operatorname{diag}(J_1, \ldots, J_m)$.
    For each Jordan block $J_k(\lambda)$:
    $$f(J_k(\lambda)) = \begin{pmatrix} f(\lambda) & f'(\lambda) & \frac{f''(\lambda)}{2!} & \cdots & \frac{f^{(k-1)}(\lambda)}{(k-1)!} \\ 0 & f(\lambda) & f'(\lambda) & \cdots & \vdots \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ 0 & 0 & \cdots & f(\lambda) & f'(\lambda) \\ 0 & 0 & \cdots & 0 & f(\lambda) \end{pmatrix}$$
    Then $f(A) = P \operatorname{diag}(f(J_1), \ldots, f(J_m)) P^{-1}$.

!!! definition "Definition 13.2 (Definition via Power Series)"
    If a scalar function $f(z)$ has a Taylor series $\sum_{k=0}^\infty c_k z^k$, and its radius of convergence is greater than the magnitudes of all eigenvalues of $A$, then:
    $$f(A) = \sum_{k=0}^\infty c_k A^k$$

---

## 13.2 The Matrix Exponential $e^A$

!!! theorem "Theorem 13.1 (Properties of the Matrix Exponential)"
    1.  **Definition**: $e^A = I + A + \frac{A^2}{2!} + \cdots$. This series converges absolutely for every square matrix $A$.
    2.  **Differentiation**: $\frac{d}{dt} e^{At} = A e^{At}$. This is the key to solving $\mathbf{x}' = A\mathbf{x}$.
    3.  **Multiplication**: If $AB = BA$, then $e^{A+B} = e^A e^B$.
    4.  **Determinant**: $\det(e^A) = e^{\operatorname{tr}(A)}$.

---

## 13.3 Computation: Interpolation Method

!!! technique "Technique: Sylvester-Lagrange Interpolation"
    If the eigenvalues of $A$ are $\lambda_1, \ldots, \lambda_k$ with maximum Jordan block sizes $n_1, \ldots, n_k$, find a polynomial $q(\lambda)$ such that its values and derivatives at $\lambda_i$ match $f(\lambda)$ and its derivatives. Then:
    $$f(A) = q(A)$$
    This avoids complex Jordan decompositions and requires only eigenvalues and powers of $A$.

---

## Exercises

1. **[Calculation] Find $e^A$ for $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$.**

   ??? success "Solution"
       $A^2 = O$, so $e^A = I + A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$.

2. **[Identity] Prove $f(P A P^{-1}) = P f(A) P^{-1}$.**

   ??? success "Solution"
       Using the power series definition: $(P A P^{-1})^k = P A^k P^{-1}$. Substituting this into the summation and factoring out $P$ and $P^{-1}$ yields the result.

3. **[Diagonal] If $A = \operatorname{diag}(\lambda_1, \lambda_2)$, find $\sin A$.**

   ??? success "Solution"
       $\sin A = \operatorname{diag}(\sin \lambda_1, \sin \lambda_2)$.

4. **[Commutativity] Prove $A f(A) = f(A) A$.**

   ??? success "Solution"
       Since $A$ commutes with all its powers $A^k$, and the series converges, it must commute with the sum of the series.

5. **[Determinant] If $\operatorname{tr}(A) = 0$, find $\det(e^A)$.**

   ??? success "Solution"
       $\det(e^A) = e^{\operatorname{tr}(A)} = e^0 = 1$.

6. **[Logarithm] For $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$, find $\log A$.**

   ??? success "Solution"
       Let $A = I + N$, where $N = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$. Since $N^2=0$, the series $\log(I+N) = N - N^2/2 + \cdots$ simplifies to $N$.

7. **[Interpolation] Use interpolation to find $f(A)$ if the only eigenvalue of $A$ is 2 with algebraic multiplicity 2.**

   ??? success "Solution"
       Find $q(\lambda) = a\lambda + b$ such that $q(2)=f(2)$ and $q'(2)=f'(2)$. Once $a, b$ are found, $f(A) = aA + bI$.

8. **[ODE] What is the solution to $\mathbf{x}' = A\mathbf{x}$ with $\mathbf{x}(0) = \mathbf{x}_0$?**

   ??? success "Solution"
       $\mathbf{x}(t) = e^{At} \mathbf{x}_0$.

9. **[Trig] Prove $\cos^2 A + \sin^2 A = I$.**

   ??? success "Solution"
       This can be verified by expanding the exponential forms: $\cos A = \frac{e^{iA}+e^{-iA}}{2}$ and $\sin A = \frac{e^{iA}-e^{-iA}}{2i}$.

10. **[Derivative] Why do derivatives appear in the function values for defective matrices?**

   ??? success "Solution"
        This reflects the "coupling" effect of matrix action. In a Jordan block, the off-diagonal 1s cause an accumulation of higher-order infinitesimals when the function is expanded, manifesting mathematically as Taylor derivative terms.

## Chapter Summary

Matrix functions are the analytic continuation of operator theory:

1.  **Consistency Principle**: Definitions via series, Jordan forms, or interpolation are equivalent within the domain of convergence, ensuring logical unity between algebra and analysis.
2.  **Exponential Core**: The matrix exponential $e^A$ is the most fundamental function, transforming linear differential evolution into pure matrix multiplication—the ultimate key to time-evolution problems.
3.  **Computational Versatility**: Interpolation and Jordan decomposition provide two complementary perspectives—the former focusing on local spectral analysis and the latter on global structural deconstruction.
