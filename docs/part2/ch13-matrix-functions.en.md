# Chapter 13: Matrix Functions

<div class="context-flow" markdown>

**Prerequisites**: Jordan Canonical Form (Ch12) · Matrix Analysis Basics (Ch14) · Power Series in Calculus

**Chapter Outline**: From Scalar to Matrix Functions → Taylor Series Definition → General Definition via Jordan Canonical Form → Sylvester-Lagrange Interpolation Method → Core Functions: Matrix Exponential ($e^A$) and Logarithm ($\ln A$) → Properties and Identities of $e^A$ → Matrix Trigonometric Functions → Applications: Analytic Solutions to Linear ODEs, Quantum Evolution Operators, and State-Transition Matrices in Control Theory

**Extension**: Matrix functions elevate arithmetic algebra to analytical algebra; they are the ultimate link between discrete matrix structure and continuous dynamical evolution—the central mathematical language for any linear system evolving over time.

</div>

Matrix functions study how to apply classical scalar functions (such as $e^x, \sin x, \ln x$) to matrix variables. This is not a simple entry-wise operation but an operator computation that preserves the algebraic structure of the matrix. Matrix functions are core tools for handling time evolution, signal propagation, and complex system responses. This chapter establishes three equivalent paths for defining matrix functions and explores the most important operator: the matrix exponential.

---

## 13.1 Methods of Definition: From Series to Jordan Form

!!! definition "Definition 13.1 (Definition via Taylor Series)"
    If a scalar function $f(z)$ has a Taylor series $\sum_{k=0}^\infty c_k z^k$ with a radius of convergence larger than the magnitudes of all eigenvalues of $A$, then:
    $$f(A) = \sum_{k=0}^\infty c_k A^k$$

!!! definition "Definition 13.2 (Definition via Jordan Form)"
    If $A = P J P^{-1}$, then $f(A) = P f(J) P^{-1}$.
    For each Jordan block $J_k(\lambda)$, the function value is defined as an upper triangular matrix composed of $f$ and its derivatives:
    $$f(J_k(\lambda)) = \begin{pmatrix} f(\lambda) & f'(\lambda) & \frac{f''(\lambda)}{2!} & \cdots & \frac{f^{(k-1)}(\lambda)}{(k-1)!} \\ 0 & f(\lambda) & f'(\lambda) & \cdots & \vdots \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ 0 & 0 & \cdots & f(\lambda) & f'(\lambda) \\ 0 & 0 & \cdots & 0 & f(\lambda) \end{pmatrix}$$

---

## 13.2 The Matrix Exponential $e^A$

!!! theorem "Theorem 13.1 (Properties of the Matrix Exponential)"
    1.  **Convergence**: For any square matrix $A$, the series $\sum \frac{A^k}{k!}$ converges absolutely.
    2.  **Differentiation**: $\frac{d}{dt} e^{At} = A e^{At}$.
    3.  **Multiplication**: If $AB = BA$, then $e^{A+B} = e^A e^B$.
    4.  **Determinant Identity**: $\det(e^A) = e^{\operatorname{tr}(A)}$ (Jacobi’s Identity).

---

## 13.3 Computation: Interpolation Method

!!! technique "Technique: Sylvester-Lagrange Interpolation"
    Let the eigenvalues of $A$ be $\lambda_1, \ldots, \lambda_m$ with multiplicities $n_1, \ldots, n_m$ in the minimal polynomial.
    There exists a unique polynomial $p(z)$ of degree less than $\sum n_i$ that matches $f(z)$ and its derivatives at each eigenvalue.
    Then: $f(A) = p(A)$.
    **Benefit**: This avoids full Jordan decomposition, requiring only eigenvalues and matrix powers.

---

## Exercises

**1. [Calculation] Find $e^A$ for $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. Observe properties: This is a nilpotent matrix, $A^2 = O$.
    2. Use the series definition: $e^A = I + A + \frac{A^2}{2!} + \cdots$.
    3. Substitute: $e^A = I + A + O = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} + \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$.
    **Conclusion**: $e^A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$.

**2. [Diagonal] Prove that if $A = \operatorname{diag}(d_1, \ldots, d_n)$, then $f(A) = \operatorname{diag}(f(d_1), \ldots, f(d_n))$.**

??? success "Solution"
    **Proof:**
    1. From the series definition: $A^k = \operatorname{diag}(d_1^k, \ldots, d_n^k)$.
    2. $f(A) = \sum c_k A^k = \sum c_k \operatorname{diag}(d_1^k, \ldots, d_n^k)$.
    3. Grouping components: $= \operatorname{diag}(\sum c_k d_1^k, \ldots, \sum c_k d_n^k)$.
    4. Each diagonal term is the series definition of $f(d_i)$.
    **Conclusion**: For diagonal matrices, the function applies to each diagonal element independently.

**3. [Commutativity] Prove that $A$ commutes with its function value: $A f(A) = f(A) A$.**

??? success "Solution"
    **Proof:**
    1. Using the series representation: $f(A) = \sum c_k A^k$.
    2. $A f(A) = A (\sum c_k A^k) = \sum c_k A^{k+1}$.
    3. $f(A) A = (\sum c_k A^k) A = \sum c_k A^{k+1}$.
    **Conclusion**: A matrix always commutes with any function of itself. This reflects that an operator commutes with the algebra it generates.

**4. [Determinant] If $\operatorname{tr}(A) = 0$, find $\det(e^A)$.**

??? success "Solution"
    **Application:**
    Using the identity $\det(e^A) = e^{\operatorname{tr}(A)}$.
    **Calculation:**
    $\det(e^A) = e^0 = 1$.
    **Physical Meaning**: This implies that if the generator of an evolution is trace-less (like a rotation), the evolution is volume-preserving.

**5. [Logarithm] Find $\ln A$ for $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$.**

??? success "Solution"
    **Method:**
    1. Let $A = I + N$, where $N = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$.
    2. Use the Mercator series: $\ln(I+N) = N - \frac{N^2}{2} + \frac{N^3}{3} - \cdots$.
    3. Since $N^2 = O$, the series truncates.
    **Conclusion**: $\ln A = N = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$.

**6. [Interpolation] Use interpolation to find $f(A)$ if $A$ has only one eigenvalue 2 with algebraic multiplicity 2.**

??? success "Solution"
    **Steps:**
    1. The minimal polynomial is $(z-2)^2$.
    2. Let the interpolation polynomial be $p(z) = a z + b$.
    3. Matching conditions: $p(2) = f(2)$ and $p'(2) = f'(2)$.
    4. Solve: $2a + b = f(2)$ and $a = f'(2)$.
    5. Result: $b = f(2) - 2f'(2)$.
    **Conclusion**: $f(A) = f'(2) A + [f(2) - 2f'(2)] I$.

**7. [ODEs] If a system satisfies $\mathbf{x}' = A\mathbf{x}$ with $\mathbf{x}(0) = \mathbf{x}_0$, what is the solution?**

??? success "Solution"
    **Conclusion:**
    $\mathbf{x}(t) = e^{At} \mathbf{x}_0$.
    **Reasoning**: Differentiating gives $\frac{d}{dt}(e^{At} \mathbf{x}_0) = A e^{At} \mathbf{x}_0 = A \mathbf{x}(t)$. At $t=0$, $e^O = I$, so it satisfies the initial condition.

**8. [Powers] Prove $(e^A)^k = e^{kA}$ for any integer $k$.**

??? success "Solution"
    **Proof:**
    Since $A$ commutes with itself, the multiplicative property gives $e^A e^A = e^{A+A} = e^{2A}$. By induction, this holds for any $k \ge 0$. For negative integers, $e^A e^{-A} = e^O = I$, so $e^{-A} = (e^A)^{-1}$.

**9. [Trig] If $A^2 = -I$, find $\cos A$.**

??? success "Solution"
    **Series expansion:**
    $\cos A = I - \frac{A^2}{2!} + \frac{A^4}{4!} - \cdots$.
    Substitute $A^2 = -I, A^4 = I \ldots$:
    $\cos A = I(1 - (-1)/2! + 1/4! - \cdots) = I(1 + 1/2! + 1/4! + \cdots)$.
    The sum is the Taylor series for $\cosh(1)$.
    **Conclusion**: $\cos A = (\cosh 1) I$.

**10. [Numerical] Why is $e^A$ difficult to compute for large matrices?**

??? success "Solution"
    **Challenges:**
    1. **Rounding Error**: Taylor series terms can have alternating signs, leading to catastrophic cancellation and loss of precision.
    2. **Computational Cost**: Standard algorithms (Scaling and Squaring) require multiple matrix multiplications, which are $O(n^3)$ each.
    3. **Conditioning**: If $A$ is close to being defective (non-normal), the computation becomes extremely sensitive to input noise.

## Chapter Summary

Matrix functions are the analytic extension of operator theory:

1.  **Consistency Principle**: Definitions via series, Jordan form, or interpolation are equivalent as long as the function is analytic on the spectrum, ensuring logical rigor.
2.  **Exponential Core**: The matrix exponential $e^A$ is the DNA of linear dynamical systems, compressing complex differential evolution into pure algebraic operator actions.
3.  **Computational Versatility**: Interpolation provides a shortcut bypassing Jordan decomposition, revealing that matrix properties are entirely determined by their local analytic behavior at eigenvalues.
