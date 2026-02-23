# Chapter 47B: Fréchet Derivatives and Higher-order Theory

<div class="context-flow" markdown>

**Prerequisites**: Matrix Calculus Basics (Ch47A) · Matrix Analysis (Ch14) · Introduction to Functional Analysis

**Chapter Outline**: From Finite-dimensional to Operator Calculus → Definition of Gateaux Derivative (Directional) → Definition of Fréchet Derivative (Total Differential) and its Linear Operator Representation → Fréchet Derivatives of Matrix Functions $L_f(A, E)$ → Key Identity: Vectorization via Kronecker Products → Higher-order Chain Rules → Derivatives of the Inverse and Matrix Exponential → Applications: Condition Number Analysis, Newton’s Method for Nonlinear Matrix Equations, and Optimization on Matrix Manifolds

**Extension**: The Fréchet derivative is the ultimate language for studying the "sensitivity" of matrix functions; it not only describes how a function value changes with input but also reveals the microscopic curvature within the operator space through its linear properties—the cornerstone of advanced numerical stability theory.

</div>

In elementary matrix calculus, we deal with gradients of scalar functions. However, when the mapping itself is matrix-to-matrix (such as $f(A) = e^A$ or $f(A) = A^{-1}$), the derivative is no longer a simple matrix but a **linear operator**. The **Fréchet derivative** is the rigorous expression of this infinitesimal linear approximation. This chapter introduces how to characterize these operators and compute them numerically using Kronecker product techniques.

---

## 47B.1 Fréchet and Directional Derivatives

!!! definition "Definition 47B.1 (Fréchet Derivative)"
    The Fréchet derivative of a mapping $f: M_n \to M_n$ at $A$ is a linear operator $L_f(A, \cdot)$ such that:
    $$f(A + E) = f(A) + L_f(A, E) + o(\|E\|)$$
    where $E$ is an infinitesimal perturbation.

!!! definition "Definition 47B.2 (Gateaux Derivative)"
    If the total differential exists, the derivative in direction $E$ can be calculated via the limit:
    $$L_f(A, E) = \lim_{h \to 0} \frac{f(A + hE) - f(A)}{h} = \left. \frac{d}{dt} f(A + tE) \right|_{t=0}$$

---

## 47B.2 Kronecker Product Representation

!!! theorem "Theorem 47B.1 (Vectorized Form of the Derivative)"
    As a linear operator, the action of the Fréchet derivative can be represented by a Kronecker product matrix $K_f(A)$ such that:
    $$\operatorname{vec}(L_f(A, E)) = K_f(A) \operatorname{vec}(E)$$
    **Significance**: This formula transforms abstract operator actions into standard matrix-vector multiplications, which is fundamental for computing matrix derivatives in numerical software.

---

## 47B.3 Fréchet Derivatives of Typical Functions

!!! example "Example 47B.1 (Inverse and Exponential)"
    1.  **Inverse** $f(A) = A^{-1}$: $L_f(A, E) = -A^{-1} E A^{-1}$.
    2.  **Square** $f(A) = A^2$: $L_f(A, E) = AE + EA$.
    3.  **Exponential** $f(A) = e^A$: The derivative is given by the integral formula $\int_0^1 e^{A(1-s)} E e^{As} ds$.

---

## Exercises

**1. [Basics] Find the Fréchet derivative of $f(A) = A^2$ using the limit definition.**

??? success "Solution"
    **Steps:**
    1. Compute $f(A+hE) = (A+hE)(A+hE) = A^2 + h(AE + EA) + h^2 E^2$.
    2. $f(A+hE) - f(A) = h(AE + EA) + h^2 E^2$.
    3. Divide by $h$ and take the limit $h \to 0$:
    **Conclusion**: $L_f(A, E) = AE + EA$. Note: Because matrices do not commute, the result is not $2AE$.

**2. [Vectorization] Express $L_f(A, E) = AE + EA$ in the Kronecker product form $K_f(A)$.**

??? success "Solution"
    **Derivation:**
    1. $\operatorname{vec}(AE) = (I \otimes A) \operatorname{vec}(E)$.
    2. $\operatorname{vec}(EA) = (A^T \otimes I) \operatorname{vec}(E)$.
    3. $\operatorname{vec}(L) = (I \otimes A + A^T \otimes I) \operatorname{vec}(E)$.
    **Conclusion**: $K_f(A) = I \otimes A + A^T \otimes I$. This is the Kronecker sum of $A$ and $A^T$.

**3. [Inverse] Prove that for $f(A) = A^{-1}$, the derivative is $L_f(A, E) = -A^{-1} E A^{-1}$.**

??? success "Solution"
    **Proof:**
    1. Start with $(A+E)(A+E)^{-1} = I$.
    2. Expand: $(A+E)(A^{-1} + L_f + \cdots) = I + E A^{-1} + A L_f + \cdots = I$.
    3. Neglecting higher-order terms, set the first-order sum to zero: $E A^{-1} + A L_f = O$.
    4. Solve for $L_f$: $A L_f = -E A^{-1} \implies L_f = -A^{-1} E A^{-1}$.

**4. [Second Derivative] Find the second-order Fréchet derivative $D^2 f(A)(E, H)$ for $f(A) = A^2$.**

??? success "Solution"
    **Derivation:**
    1. The first derivative is $L(A, E) = AE + EA$.
    2. Differentiate $L$ with respect to $A$ in direction $H$:
    3. $\delta L = (A+H)E + E(A+H) - (AE + EA) = HE + EH$.
    **Conclusion**: $D^2 f(A)(E, H) = HE + EH$. In this case, the second derivative is independent of $A$ (as the original function was quadratic).

**5. [Determinant] Calculate the Fréchet derivative of $\det(A)$ (assume $A$ is invertible).**

??? success "Solution"
    **Jacobi’s Formula:**
    $\delta \det(A) = \det(A) \operatorname{tr}(A^{-1} E)$.
    Since the result is a scalar, this is equivalent to the gradient $\nabla \det A = \det(A) (A^{-1})^T$ from Ch47A in the inner product sense.

**6. [Chain Rule] If $h(A) = f(g(A))$, write the composition of their Fréchet derivatives.**

??? success "Solution"
    **Formula:**
    $L_h(A, E) = L_f(g(A), L_g(A, E))$.
    **Significance**: This shows that the composition of derivative operators follows the nesting of linear operators, not simple matrix multiplication.

**7. [Exponential] Why isn't the derivative of $e^A$ simply $e^A E$?**

??? success "Solution"
    **Reasoning:**
    The exponential identity $e^{A+E} = e^A e^E$ only holds if $A$ and $E$ commute. In general, they do not, so the terms in the Taylor expansion cannot be rearranged freely. One must use the integral representation (as in Example 47B.1).

**8. [Application] What is the "Condition Number" of a matrix function?**

??? success "Solution"
    **Definition:**
    $\kappa_f(A) = \frac{\|L_f(A, \cdot)\| \cdot \|A\|}{\|f(A)\|}$, where $\|L_f\|$ is the operator norm of the linear derivative.
    It measures the relative sensitivity of the matrix function computation to perturbations in the input.

**9. [Calculation] Find the derivative of $f(A) = A^3$.**

??? success "Solution"
    **Derivation:**
    $\delta(AAA) = (\delta A)AA + A(\delta A)A + AA(\delta A)$.
    **Conclusion**: $L_f(A, E) = EAA + AEA + AAE$.

**10. [Numerical] Briefly explain the Complex Step method for approximating Fréchet derivatives.**

??? success "Solution"
    **Idea:**
    Use $f(A + i h E) \approx f(A) + i h L_f(A, E)$.
    Take the imaginary part: $L_f(A, E) \approx \frac{\operatorname{Im}(f(A + i h E))}{h}$.
    This method avoids subtraction cancellation errors found in finite differences, providing near-machine precision for the derivative.

## Chapter Summary

The Fréchet derivative is the total differential tool in operator space:

1.  **Essence of Linear Approximation**: It simplifies complex non-linear matrix functions locally into linear operators acting on perturbation matrices, establishing the standard language for sensitivity analysis.
2.  **Vectorization Bridge**: The introduction of the Kronecker matrix $K_f(A)$ realizes the jump from abstract operator theory to concrete numerical calculation, serving as the basis for high-performance math libraries.
3.  **Precision in Non-commutativity**: By deriving derivatives for core functions like $A^{-1}$ and $e^A$, this chapter highlights the deepest difference between matrix and scalar calculus—the strict respect for order and non-commutativity.
