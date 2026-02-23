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

****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ## Chapter Summary

Matrix functions are the analytic continuation of operator theory:


****: Definitions via series, Jordan forms, or interpolation are equivalent within the domain of convergence, ensuring logical unity between algebra and analysis.

****: The matrix exponential $e^A$ is the most fundamental function, transforming linear differential evolution into pure matrix multiplication—the ultimate key to time-evolution problems.

****: Interpolation and Jordan decomposition provide two complementary perspectives—the former focusing on local spectral analysis and the latter on global structural deconstruction.
