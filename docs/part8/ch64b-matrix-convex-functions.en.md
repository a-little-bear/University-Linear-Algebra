# Chapter 64B: Matrix Convex Functions

<div class="context-flow" markdown>

**Prerequisites**: Matrix Analysis (Ch14) · Positive Definite Matrices (Ch16) · Convex Sets (Ch64A) · Matrix Inequalities (Ch18)

**Chapter Outline**: Definition of Matrix Convexity → Difference between Operator Convex and Operator Monotone → Key Operator Convex Functions ($x^2, x^{-1}, x\log x$) → Jensen’s Operator Inequality → Kraus’s Theorem (Representation Theory) → Basics of Löwner’s Theorem → Strong Subadditivity (SSA) & Quantum Entropy → Applications: Convexity of Matrix Means and Operator Divergence

**Extension**: Matrix convexity is a far stronger constraint than scalar convexity; a function can be convex on the real line but fail to be convex in the matrix sense. This theory reveals how non-commutative operators preserve order under non-linear mappings, forming the cutting edge of information theory and functional analysis.

</div>

In classical analysis, a convex function satisfies $f(\lambda x + (1-\lambda)y) \le \lambda f(x) + (1-\lambda)f(y)$. However, when we replace the variables with matrices, the situation becomes significantly more complex. **Matrix Convex Functions** (also known as operator convex functions) require this inequality to hold in the sense of the **Löwner partial order**. This stringent requirement excludes most ordinary convex functions, leaving only a class of operators with deep analytic backgrounds (such as Pick functions). This chapter explores these powerful non-linear operator properties.

---

## 64B.1 Definition of Operator Convexity

!!! definition "Definition 64B.1 (Operator Convex Function)"
    Let $f$ be a real-valued function on an interval $I$. $f$ is **operator convex** if for all symmetric matrices $A, B$ with eigenvalues in $I$, and any $\lambda \in [0, 1]$:
    $$f(\lambda A + (1-\lambda)B) \preceq \lambda f(A) + (1-\lambda)f(B)$$
    The inequality is interpreted in the Löwner order ($X \preceq Y \iff Y - X$ is positive semi-definite).

!!! warning "Operator Convex vs. Scalar Convex"
    Operator convexity implies scalar convexity, but the converse is false. For example, $f(x) = x^4$ is convex on $\mathbb{R}$, but it is **not** operator convex.

---

## 64B.2 Fundamental Operator Convex Functions

!!! theorem "Theorem 64B.1 (Key Examples)"
    1.  **Inverse**: $f(x) = 1/x$ is operator convex on $(0, \infty)$ (and operator monotone decreasing).
    2.  **Square**: $f(x) = x^2$ is operator convex on $\mathbb{R}$.
    3.  **Operator Entropy**: $f(x) = x \log x$ is operator convex on $(0, \infty)$ (the foundation of quantum information theory).
    4.  **Negative Powers**: $f(x) = x^r$ is operator convex for $r \in [1, 2]$ or $r \in [-1, 0]$.

---

## 64B.3 Jensen’s Operator Inequality

!!! theorem "Theorem 64B.2 ( Jensen’s Operator Inequality)"
    Let $f$ be an operator convex function, and let $\sum V_i^* V_i = I$ (a partition of unity). Then:
    $$f\left( \sum V_i^* A_i V_i \right) \preceq \sum V_i^* f(A_i) V_i$$
    **Application**: This is the core tool for proving data processing inequalities for quantum channels.

---

## 64B.4 Kraus’s Theorem

!!! theorem "Theorem 64B.3 (Kraus’s Theorem)"
    A function $f$ is operator convex if and only if it can be represented in the following integral form:
    $$f(x) = a + bx + cx^2 + \int \frac{(x-\lambda)^2}{1+\lambda^2} d\mu(\lambda)$$
    where $c \ge 0$ and $\mu$ is a positive measure. This links operator convexity to the theory of analytic continuation in complex analysis.

---

## Exercises

1.  **[Basics] Prove that $f(x) = ax+b$ is both operator convex and operator concave.**
    ??? success "Solution"
        Since $f(\lambda A + (1-\lambda)B) = a(\lambda A + (1-\lambda)B) + bI = \lambda f(A) + (1-\lambda)f(B)$, the equality holds identically.

2.  **[Counter-example] Is $f(x) = x^3$ operator convex on the set of positive definite matrices?**
    ??? success "Solution"
        No. While it is scalar convex on $\mathbb{R}^+$, it fails the operator convexity condition for matrices of order $n \ge 2$.

3.  **[Property] Prove: If $f$ is operator convex, then $g(x) = f(x) + c$ is also operator convex.**
    ??? success "Solution"
        The constant terms cancel out on both sides of the inequality, leaving the Löwner relationship unchanged.

4.  **[Inverse] Explain why $A \preceq B \implies B^{-1} \preceq A^{-1}$ for positive definite matrices.**
    ??? success "Solution"
        This follows from the fact that $f(x) = 1/x$ is an operator monotone decreasing function.

5.  **[Sum] Is the sum of two operator convex functions necessarily operator convex?**
    ??? success "Solution"
        Yes. The Löwner order is closed under addition.

6.  **[Entropy] In quantum mechanics, why is the von Neumann entropy $S(\rho) = -\operatorname{tr}(\rho \log \rho)$ concave?**
    ??? success "Solution"
        Because $f(x) = x \log x$ is operator convex, so its negative is operator concave. The trace operator preserves this concavity.

7.  **[Monotonicity] Is an operator convex function always operator monotone?**
    ??? success "Solution"
        No. For example, $x^2$ is operator convex but is not operator monotone on $(0, \infty)$ (unless restricted to a specific domain).

8.  **[Calculation] Let $A = \operatorname{diag}(1, 2)$ and $B = \operatorname{diag}(2, 1)$. Verify operator convexity for $f(x)=x^2$.**
    ??? success "Solution"
        $(\frac{A+B}{2})^2 = (1.5 I)^2 = 2.25 I$. $\frac{A^2+B^2}{2} = \frac{\operatorname{diag}(1, 4) + \operatorname{diag}(4, 1)}{2} = 2.5 I$. Since $2.25 I \preceq 2.5 I$, convexity holds.

9.  **[SSA] What is Strong Subadditivity (SSA) of quantum entropy?**
    ??? success "Solution"
        $S(\rho_{ABC}) + S(\rho_B) \le S(\rho_{AB}) + S(\rho_{BC})$. Its algebraic proof relies heavily on the operator convexity of $x \log x$ and the resulting trace inequalities.

10. **[Differentiation] Why does operator convexity require $f$ to be at least $C^2$?**

   ??? success "Solution"
        Actually, operator convexity is a much stronger condition; it implies that the second Fréchet derivative, as a bilinear form, must be positive definite (or that the second-order divided difference matrix is PSD).

## Chapter Summary

Matrix convex functions establish stability criteria under non-linear operator actions:

1.  **Rigidity of Order**: Operator convexity proves that in non-commutative environments, non-linear transformations that preserve partial ordering are extremely rare, setting strict mathematical boundaries for energy evolution in physical systems.
2.  **Analytic Link**: Theorems by Kraus and Löwner link purely algebraic inequalities to analytic properties in the complex plane, revealing that operator convexity is essentially a projection of certain holomorphic structures.
3.  **Bedrock of Information**: Operator convex functions, epitomized by $x \log x$, support the entirety of entropy principles and stability analysis in quantum information theory, proving that algebraic convexity is the fundamental logic of cosmic information conservation.
