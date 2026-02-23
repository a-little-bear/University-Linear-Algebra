# Chapter 64B: Matrix Convex Functions and Operator Monotonicity

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Matrix Inequalities (Ch18) · Convex Sets (Ch64A) · Majorization (Ch31)

**Chapter Outline**: Scalar Convex Matrix Functions → Matrix Concave Functions → Operator Convex Functions (Matrix Order) → Löwner's Theorem (Operator Monotone $\leftrightarrow$ Pick Functions) → Hansen-Pedersen Theorem → Lieb's Concavity Theorem → Matrix Perspective Functions → Schur Convexity

**Extension**: Operator convex/monotone theory connects matrix analysis to functional analysis and quantum information (relative entropy, data processing inequalities).

</div>

The theory of matrix-valued functions involves two distinct levels: scalar convexity (where the function value is a scalar) and the deeper operator convexity (where inequalities hold in the Löwner partial order). Operator convexity is a significantly more restrictive condition, revealing the deep topological and analytical structure of the space of operators.

---

## 64B.1 Two Dimensions of Convexity

!!! definition "Definition 64B.1 (Scalar vs. Operator Convexity)"
    - **Scalar Convex**: $f(tA + (1-t)B) \le t f(A) + (1-t) f(B)$ (functional value is a scalar).
    - **Operator Convex**: $f(tA + (1-t)B) \preceq t f(A) + (1-t) f(B)$ (inequality holds in the Löwner order).

!!! theorem "Theorem 64B.5 (Löwner's Theorem, 1934)"
    A function $f$ is operator monotone on an interval iff it can be analytically extended to a Pick function (mapping the upper half-plane to itself).

---

## Exercises

1. **[Fundamentals] Prove that $f(A) = \lambda_{\max}(A)$ is a (scalar) convex function on symmetric matrices.**
   ??? success "Solution"
       $\lambda_{\max}(A) = \max_{\|x\|=1} x^T A x$. Since for each $x$, the map $A \mapsto x^T Ax$ is linear (and thus convex), $\lambda_{\max}$ is convex as the pointwise supremum of a family of convex functions.

2. **[Calculation] Is $f(t) = t^2$ operator convex on $[0, \infty)$? Is it operator monotone?**
   ??? success "Solution"
       - **Operator Convex**: Yes, the square of a matrix is operator convex.
       - **Operator Monotone**: No. Counterexample: $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} \preceq \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix} = B$, but $A^2 \npreceq B^2$ (as $\det(B^2 - A^2) = -4 < 0$).

3. **[Concavity] Show that the log-determinant $f(A) = \log \det A$ is a (scalar) concave function on $S_n^{++}$.**
   ??? success "Solution"
       Using the eigenvalue representation, $\log \det (\frac{A+B}{2}) = \sum \log \lambda_i(\frac{A+B}{2})$. The concavity of the $\log$ function combined with the Minkowski inequality for determinants proves scalar concavity.

4. **[Ky Fan] Define the Ky Fan $k$-norm and explain why it is a convex function.**
   ??? success "Solution"
       The Ky Fan $k$-norm is the sum of the $k$ largest singular values. It can be characterized as the maximum of $\operatorname{tr}(U^T AV)$ over the Stiefel manifold ($U^T U = V^T V = I_k$). As the supremum of a family of linear functions, it is convex.

5. **[Hansen-Pedersen] State the operator version of Jensen's Inequality.**
   ??? success "Solution"
       For an operator convex function $f$, $f(C^* AC) \preceq C^* f(A) C$ holds for all matrices $A$ with spectrum in the domain and all contractions $C$ ($C^* C \preceq I$).

6. **[Lieb's Theorem] What is the core application of Lieb's Concavity Theorem in quantum information?**
   ??? success "Solution"
       It is the fundamental tool used to prove the **Strong Subadditivity** of von Neumann entropy, ensuring the physical consistency of information flow in quantum systems.

7. **[Matrix Means] Prove that the matrix geometric mean $A \# B$ is jointly concave in $(A, B)$.**
   ??? success "Solution"
       Using the characterization $A \# B = \max \{ X : \begin{pmatrix} A & X \\ X & B \end{pmatrix} \succeq 0 \}$, the joint concavity follows because the feasible set of the LMI is convex in the block parameters $(A, B)$.

8. **[Schur Convexity] Determine the Schur convexity of $\phi(x) = \sum x_i^2$.**
   ??? success "Solution"
       Since $\phi$ is symmetric and the component function $t^2$ is convex, $\phi$ is Schur convex. This implies that for a fixed trace, the sum of squares is maximized when the eigenvalues are most "spread out."

9. **[Perspective] If $f(x) = 1/x$, is its perspective function $g(A, t) = t(A/t)^{-1} = t^2 A^{-1}$ jointly convex in $(A, t)$?**
   ??? success "Solution"
       Yes. Since $1/x$ is an operator convex function, its perspective $t f(A/t)$ is jointly convex in the Löwner order on $S_n^{++} \times (0, \infty)$.

10. **[Trace Functions] Contrast the convexity of $\operatorname{tr}(f(A))$ with the operator convexity of $f(A)$.**
    ??? success "Solution"
        If $f$ is a convex scalar function, then $\operatorname{tr}(f(A))$ is always a convex scalar function of $A$. However, operator convexity of $f(A)$ requires much stricter conditions (analytic Pick functions). Every operator convex function induces a convex trace function, but not vice versa (e.g., $t^3$ on $[0, \infty)$).

## Chapter Summary

This chapter details the analytical and ordering properties of matrix-valued functions:

1. **Hierarchy of Convexity**: Distinguished between scalar convexity and the much more rigid operator convexity.
2. **Analytical Extensions**: Explored Löwner's theorem, linking operator monotonicity to the theory of Pick functions in complex analysis.
3. **Operator Jensen Inequalities**: Extended classical inequalities to the operator domain via the Hansen-Pedersen framework.
4. **Information Metrics**: Analyzed Lieb's concavity and matrix perspectives, establishing the mathematical foundation for quantum entropy and information geometry.
