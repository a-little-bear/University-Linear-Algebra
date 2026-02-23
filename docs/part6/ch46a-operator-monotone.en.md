# Chapter 46A: Operator Monotone Functions

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Matrix Analysis (Ch14) · Löwner Partial Order (Ch16)

**Chapter Outline**: From Scalar to Operator Monotonicity → Definition of Operator Monotone Functions → Operator Convex Functions → Core Theorem: Loewner’s Theorem (Analytic Continuation and Pick Functions) → Classic Examples (Square Root, Logarithm) → Counterexamples: Why $t^2$ is not Operator Monotone → Transitivity of Operator Inequalities → Applications: Quantum Information (Concavity of Entropy), Foundation of Matrix Means, and Stability Estimates in Statistical Physics

**Extension**: Operator monotonicity studies the preservation of matrix partial orders under non-linear mappings; it is one of the most advanced topics in matrix analysis, revealing the incredible unity between complex analytic continuation and the algebraic order structure of matrices. It is the algebraic core for constructing Matrix Means (Ch46B).

</div>

In scalar analysis, if $f(t)$ is monotonically increasing and $a \ge b$, then $f(a) \ge f(b)$. In the matrix world, however, this intuition often fails. Even if $A \succeq B \succeq 0$, it does not necessarily follow that $A^2 \succeq B^2$. **Operator Monotone Functions** are those specific non-linear functions that preserve the **Löwner partial order** of matrices. This chapter introduces the elegant theory established by Loewner, revealing how the analyticity of a function dictates matrix inequalities.

---

## 46A.1 Definitions and Core Classification

!!! definition "Definition 46A.1 (Operator Monotone Function)"
    A function $f$ defined on an interval $I$ is **operator monotone** if for any Hermitian matrices $A, B$ with eigenvalues in $I$:
    $$A \succeq B \Rightarrow f(A) \succeq f(B)$$

!!! definition "Definition 46A.2 (Operator Convex Function)"
    A function $f$ is **operator convex** if:
    $$f(\lambda A + (1-\lambda)B) \preceq \lambda f(A) + (1-\lambda)f(B)$$
    for all $A, B$ and $\lambda \in [0, 1]$. This requires the function's action on a mixture of operators to stay "below" the mixture of the function values.

---

## 46A.2 Core Theorem: Loewner's Theorem

!!! theorem "Theorem 46A.1 (Loewner's Theorem)"
    A function $f$ is operator monotone on $(a, b)$ iff it can be analytically continued to the complex upper half-plane and maps the upper half-plane into itself (a Pick function).
    **Significance**: This binds an algebraic property (order preservation) to a geometric analytic property (region mapping).

---

## 46A.3 Classic Examples and Counterexamples

!!! example "Example 46A.1 (Classic Operator Monotone Functions)"
    1.  **Power Functions**: $f(t) = t^p$ is operator monotone for $0 \le p \le 1$ (e.g., the square root $\sqrt{t}$).
    2.  **Logarithm**: $f(t) = \ln t$ is operator monotone on $(0, \infty)$.
    3.  **Linear Fractional**: $f(t) = -1/t$ is operator monotone on $(0, \infty)$.

!!! warning "Counterexample: The Failure of $t^2$"
    Even if $A \succeq B \succeq 0$, it is generally **not** true that $A^2 \succeq B^2$. Squaring operations can "warp" the original order structure of the matrix space.

---

## Exercises

**1. [Basics] If $A = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$ and $B = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$, verify if $A^2 \succeq B^2$.**

??? success "Solution"
    **Calculation:**
    1. $A - B = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} \succeq 0$, so $A \succeq B$.
    2. $A^2 = \begin{pmatrix} 4 & 0 \\ 0 & 1 \end{pmatrix}, B^2 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$.
    3. $A^2 - B^2 = \begin{pmatrix} 3 & 0 \\ 0 & 0 \end{pmatrix} \succeq 0$.
    **Conclusion**: It holds in this diagonal case. However, for non-diagonal matrices, this is not a general law.

**2. [Counterexample] Provide an example where $A \succeq B \ge 0$ but $A^2 \nsucceq B^2$.**

??? success "Solution"
    **Construction:**
    Let $B = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$ and $A = B + \begin{pmatrix} \epsilon & 0 \\ 0 & 0 \end{pmatrix}$.
    Under the influence of cross-terms, $A^2 - B^2$ can develop negative eigenvalues. This proves $f(t)=t^2$ is not monotone in the operator sense.

**3. [Properties] Prove: If $f$ is operator monotone, it must be scalar monotone.**

??? success "Solution"
    **Proof:**
    1. Consider $1 \times 1$ matrices (scalars) $a, b$.
    2. The definition of operator monotonicity requires $a \ge b \Rightarrow f(a) \ge f(b)$.
    3. This is exactly the definition of scalar monotonicity.
    **Conclusion**: Operator monotonicity is a much stricter constraint than scalar monotonicity.

**4. [Square Root] Prove: If $A \succeq B \succeq 0$, then $\sqrt{A} \succeq \sqrt{B}$.**

??? success "Solution"
    **Reasoning:**
    The function $f(t) = t^{1/2}$ satisfies Loewner’s theorem on $[0, \infty)$ (its analytic branch maps the upper half-plane to the upper half-plane). Thus, it is operator monotone and preserves the matrix order.

**5. [Inverse] Prove $f(t) = -1/t$ is operator monotone on $(0, \infty)$.**

??? success "Solution"
    **Proof:**
    1. Let $A \succeq B \succ 0$.
    2. Matrix properties (Ch16) show inversion reverses the order: $B^{-1} \succeq A^{-1}$.
    3. Multiplying by $-1$ reverses it again: $-A^{-1} \succeq -B^{-1}$.
    4. Thus $f(A) \succeq f(B)$.

**6. [Convexity] Is $f(t) = t^2$ operator convex?**

??? success "Solution"
    **Yes.**
    **Proof Sketch**: Expanding $(\lambda A + (1-\lambda)B)^2$ and utilizing $A^2, B^2$, the middle term $AB+BA$ is always controlled by a combination of $A^2$ and $B^2$. This shows that while $t^2$ does not preserve order, it does preserve "convexity."

**7. [Application] Why are entropy properties in quantum mechanics related to operator monotonicity?**

??? success "Solution"
    Quantum relative entropy involves $f(t) = \ln t$ or $f(t) = t \ln t$. Since $\ln t$ is operator monotone (and $t \ln t$ is operator convex), these properties guarantee that information content does not spontaneously increase under mixing, ensuring consistency with the Second Law of Thermodynamics.

**8. [Calculation] Compute $f(A)$ where $f(t) = \sqrt{t}$ and $A = \begin{pmatrix} 5 & 4 \\ 4 & 5 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. Eigen-decomposition: Eigenvalues are 9 (with vector $(1,1)$) and 1 (with vector $(1,-1)$).
    2. Take roots of eigenvalues: 3 and 1.
    3. Reconstruct: $\sqrt{A} = \frac{1}{2} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix} \begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$.

**9. [Composition] If $f$ and $g$ are operator monotone, is $f \circ g$ also operator monotone?**

??? success "Solution"
    **Yes.**
    If $A \succeq B$, then $g(A) \succeq g(B)$ by monotonicity of $g$. Then $f(g(A)) \succeq f(g(B))$ by monotonicity of $f$. The property holds.

**10. [Limits] As $p \to 0$, what does the function $f(t) = (t^p - 1)/p$ approach?**

??? success "Solution"
    **Conclusion:**
    It approaches $\ln t$.
    Since $(t^p-1)/p$ is operator monotone for each $p \in (0, 1)$, the limit $\ln t$ inherits this property. This demonstrates the closure of the class of operator monotone functions.

## Chapter Summary

Operator monotone functions are the strict "order keepers" of matrix analysis:

1.  **Guards of Partial Order**: They define which non-linear transformations respect the Loewner order of matrix space, filtering out functions like "square" which, though simple, are destructive in higher dimensions.
2.  **Analytic Link**: Loewner's theorem reveals the deep correspondence between matrix inequalities on the real line and Pick analytic functions in the complex plane—a paradigm of cross-domain unification in mathematics.
3.  **Physical Bedrock**: As properties of core operators like logs and powers, operator monotonicity provides the ultimate algebraic guarantee for entropy properties, the validity of matrix mean definitions, and system stability.
