# Chapter 64B: Matrix Convex Functions

<div class="context-flow" markdown>

**Prerequisites**: Operator Monotone Functions (Ch46A) · Convex Sets in Matrix Spaces (Ch64A) · Matrix Analysis (Ch14)

**Chapter Outline**: From Scalar to Operator Convexity → Definition of Matrix Convex Functions → Matrix Jensen Inequality → The Deep Link between Matrix Convexity and Operator Monotonicity → Typical Matrix Convex Functions: Inverse $X^{-1}$, Negative Logarithm $-\ln X$, and Power Functions $X^p$ → Trace Convexity → Applications: Fisher Information in Statistics, Free Energy Minimization in Quantum Systems, and Complexity Estimation for SDP Algorithms

**Extension**: Matrix convexity is the "license" for optimization theory to enter the matrix domain; it proves that complex matrix mappings can maintain the downward curvature of energy surfaces, serving as the underlying logic for the absolute convergence of Newton and interior-point methods in large-scale scientific computing.

</div>

In scalar calculus, a convex function $f(\lambda x + (1-\lambda)y) \le \lambda f(x) + (1-\lambda)f(y)$ is the guarantee for finding global minima. In the matrix world, we are concerned with mappings that satisfy the convexity definition under the **Löwner partial order**. **Matrix Convex Functions** describe how operators maintain the growth of "energy" or "uncertainty" under mixed inputs. This chapter introduces the criteria for matrix convexity and its central role in information theory and optimal design.

---

## 64B.1 Definition and Matrix Jensen Inequality

!!! definition "Definition 64B.1 (Matrix Convex Function)"
    A function $f$ defined on an interval $I$ is **matrix convex** if for any Hermitian matrices $A, B$ with eigenvalues in $I$:
    $$f(\lambda A + (1-\lambda)B) \preceq \lambda f(A) + (1-\lambda)f(B)$$
    for all $\lambda \in [0, 1]$.

!!! theorem "Theorem 64B.1 (Matrix Jensen Inequality)"
    If $f$ is a matrix convex function, then for any set of isometry matrices $\{V_i\}$ such that $\sum V_i^* V_i = I$:
    $$f\left( \sum V_i^* A_i V_i \right) \preceq \sum V_i^* f(A_i) V_i$$
    This property is extremely powerful when dealing with weighted averages of operators.

---

## 64B.2 Typical Matrix Convex Functions

!!! note "Common Functions"
    1.  **Inverse**: $f(X) = X^{-1}$ is matrix convex on $X \succ 0$.
    2.  **Negative Log**: $f(X) = -\ln X$ is matrix convex on $X \succ 0$.
    3.  **Power Functions**: $f(X) = X^p$ is matrix convex if $1 \le p \le 2$ or $-1 \le p \le 0$.

---

## 64B.3 Trace Convexity

!!! technique "Technique: The Trace Advantage"
    Sometimes a function is not matrix convex, but its trace $\operatorname{tr}(f(X))$ is convex. For example, $f(X) = X^p$ is **trace convex** for any $p \ge 1$. This is used in statistical physics to bound the lower limit of free energy.

---

## Exercises

**1. [Basics] Prove that $f(X) = X^2$ is a matrix convex function.**

??? success "Solution"
    **Proof:**
    1. Compute $f(\lambda A + (1-\lambda)B) = (\lambda A + (1-\lambda)B)^2$.
    2. Expand: $\lambda^2 A^2 + (1-\lambda)^2 B^2 + \lambda(1-\lambda)(AB + BA)$.
    3. Compare with $\lambda A^2 + (1-\lambda)B^2$.
    4. The difference is: $\lambda(1-\lambda)(A^2 + B^2 - AB - BA) = \lambda(1-\lambda)(A-B)^2$.
    5. Since $(A-B)^2 \succeq 0$ and $\lambda(1-\lambda) \ge 0$.
    **Conclusion**: The difference is positive semi-definite, so $f(X) = X^2$ is matrix convex.

**2. [Comparison] Why is $f(X) = X^2$ matrix convex but not operator monotone?**

??? success "Solution"
    **Nuance:**
    - **Convexity** reflects the "curvature" of the function graph. The square function always curves upward, in both scalar and matrix spaces.
    - **Monotonicity** reflects whether the function "preserves order." In non-commutative spaces, the squaring operation introduces cross-terms $AB+BA$ that can disrupt the original $A \succeq B$ order.
    **Conclusion**: This is a key insight in matrix analysis: convexity is more easily preserved than monotonicity in operator spaces.

**3. [Calculation] Determine the convexity of $f(X) = X^{-1}$ on positive matrices.**

??? success "Solution"
    **Using Schur Complement:**
    1. Consider the block matrix $M = \begin{pmatrix} X & I \\ I & f(X) \end{pmatrix}$.
    2. If $f(X) = X^{-1}$, the Schur complement is $X^{-1} - X^{-1} = 0$, placing it at the boundary of positive definiteness.
    3. Utilizing the properties of convex combinations of block matrices, one can prove $X^{-1}$ satisfies the definition of matrix convexity.
    **Conclusion**: Inversion is a highly non-linear matrix convex operation.

**4. [Application] Briefly state the significance of matrix convexity in Fisher Information matrices.**

??? success "Solution"
    In statistical estimation, the Fisher Information matrix $I(\theta)$ measures the information content of observations. Since inversion is convex, the Cramér-Rao Lower Bound ($I^{-1}$) exhibits convexity under mixture experiments, ensuring that increasing samples always (in a convexity sense) reduces estimation error.

**5. [Trace] Determine if $\phi(X) = \operatorname{tr}(X^3)$ is convex for $X \succeq 0$.**

??? success "Solution"
    **Conclusion: Yes.**
    While $X^3$ is not necessarily convex in the operator sense, the trace eliminates the non-symmetric effects of cross-terms, making $\operatorname{tr}(X^3)$ exhibit perfect scalar convexity on the positive axis.

**6. [Log] Prove $-\ln X$ is matrix convex.**

??? success "Solution"
    **Reasoning:**
    We know $\ln X$ is operator monotone (Ch46A). According to operator theory, every positive operator monotone function on $(0, \infty)$ is operator concave. Since $\ln X$ is concave, its negative $-\ln X$ is **matrix convex**.

**7. [Jensen] Use the Matrix Jensen Inequality to prove $f(U^* A U) = U^* f(A) U$ for unitary $U$.**

??? success "Solution"
    **Proof:**
    When $\sum V_i^* V_i = I$ simplifies to a single term $U^* U = I$, the Jensen inequality holds with equality (as it is a change of basis). This proves that matrix convexity is fully compatible with unitary invariance.

**8. [Stability] Why is matrix convexity required for objective functions when solving LMIs?**

??? success "Solution"
    **Reasoning:**
    A local optimum is a global optimum only if the objective function is convex and the feasible region (defined by the LMI) is a convex set. This allows us to confidently use interior-point methods to find robust controllers in large-scale spaces.

**9. [Property] If $f$ is matrix convex and $f(0) \le 0$, does $f(AXA^*) \preceq A f(X) A^*$ always hold?**

??? success "Solution"
    No, not necessarily. This is a subtle modification of the Jensen condition. Matrix convexity imposes strict constraints on the "coefficients" (weights); they must satisfy normalization conditions (like being isometries).

**10. [Application] What is the Strong Subadditivity of Von Neumann Entropy?**

??? success "Solution"
    **Connection:**
    This is one of the most profound inequalities in quantum information theory. Its proof relies fundamentally on the matrix convexity of $f(X) = X \ln X$ and its interaction with operator monotonicity. It ensures that correlations in composite quantum systems have physically reasonable bounds.

## Chapter Summary

Matrix convex functions represent the intersection of operator analysis and optimization theory:

1.  **Guarantee of Global Optimality**: They elevate the intuition of convexity from scalar space to operator partial order space, providing theoretical shielding for complex matrix-constrained optimization.
2.  **Interplay of Monotony and Convexity**: Matrix convexity is not an isolated property but is deeply complementary to operator monotonicity via analytic continuation theorems (Loewner), revealing the deep consistency of positive operator mappings.
3.  **Measure of Information**: From inverses to logs, matrix convexity provides the unique algebraic framework for characterizing statistical information, physical energy, and quantum correlation—the foundation of modern system science.
