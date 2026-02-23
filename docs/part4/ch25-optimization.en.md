# Chapter 25: Optimization and Eigenvalue Problems

<div class="context-flow" markdown>

**Prerequisites**: Matrix Calculus (Ch47A) · Positive Definiteness (Ch16) · Eigenvalues (Ch6) · Convex Sets (Ch64A)

**Chapter Outline**: Extremum Problems for Quadratic Forms → Rayleigh Quotient → Courant-Fischer Min-Max Theorem → Ky Fan's Maximum Principle → Optimization of Trace Functions → Log-barrier and Determinants → Semidefinite Programming (SDP) → Interior Point Methods for Matrices → Procrustes Problems

**Extension**: Eigenvalue problems are the solutions to the most fundamental optimization problems in geometry (maximizing variance or minimizing energy).

</div>

Many optimization problems in science and engineering can be reduced to finding the extreme values of a matrix-valued function. The most famous example is the **Rayleigh Quotient**, whose maximum and minimum values are precisely the largest and smallest eigenvalues of the matrix. This chapter bridges the gap between linear algebra and optimization, exploring how spectral theorems provide the closed-form solutions to complex multivariable extremum problems. We also introduce **Semidefinite Programming** (SDP), which generalizes linear programming to matrix variables.

---

## 25.1 Variational Principles for Eigenvalues

!!! definition "Definition 25.1 (Rayleigh Quotient)"
    For a symmetric matrix $A$, the Rayleigh quotient is:
    $$R(x) = \frac{x^T A x}{x^T x}, \quad x \neq 0$$
    The maximum and minimum of $R(x)$ are $\lambda_{\max}(A)$ and $\lambda_{\min}(A)$ respectively.

!!! theorem "Theorem 25.1 (Courant-Fischer Theorem)"
    The $k$-th largest eigenvalue of a symmetric matrix $A$ is given by:
    $$\lambda_k = \max_{\dim S=k} \min_{x \in S, x \neq 0} \frac{x^T A x}{x^T x}$$
    This characterizes all eigenvalues as solutions to nested min-max optimization problems.

---

## Exercises

1. **[Fundamentals] Maximize $x^2 + 2y^2 + 3z^2$ subject to $x^2 + y^2 + z^2 = 1$.**
   ??? success "Solution"
       This is the maximum of the Rayleigh quotient for $A = \operatorname{diag}(1, 2, 3)$. The maximum is $\lambda_{\max} = 3$, achieved at $(0, 0, 1)^T$.

2. **[Ky Fan] State Ky Fan's maximum principle for the sum of eigenvalues.**
   ??? success "Solution"
       $\sum_{i=1}^k \lambda_i(A) = \max \operatorname{tr}(U^T AU)$ where $U^T U = I_k$. The sum of the largest $k$ eigenvalues is the maximum trace achievable by projecting into a $k$-dimensional subspace.

3. **[SDP] Write the standard form of a Semidefinite Program.**
   ??? success "Solution"
       $\min \operatorname{tr}(CX)$ subject to $\operatorname{tr}(A_i X) = b_i$ and $X \succeq 0$. This generalizes linear programming by replacing vector variables with PSD matrix variables.

4. **[Procrustes] What is the Orthogonal Procrustes problem?**
   ??? success "Solution"
       Given $A$ and $B$, find an orthogonal matrix $Q$ that minimizes $\|AQ - B\|_F$. The solution is $Q = UV^T$ from the SVD of $A^T B$.

5. **[Log-Barrier] Why is $-\log\det X$ used as a barrier function in optimization?**
   ??? success "Solution"
       As $X$ approaches the boundary of the PSD cone (i.e., becomes singular), $\det X \to 0$, so $-\log\det X \to \infty$. This penalizes points near the boundary, keeping the optimization in the interior ($X \succ 0$).

6. **[Trace Optimization] Maximize $\operatorname{tr}(AX)$ subject to $X^T X = I$.**
   ??? success "Solution"
       If $A = U \Sigma V^T$, the optimal $X = VU^T$. This is a fundamental result in matrix alignment and shape analysis.

7. **[Convexity] Prove that $\lambda_{\max}(A)$ is a convex function of $A$.**
   ??? success "Solution"
       $\lambda_{\max}(A) = \max_{\|x\|=1} x^T Ax$. Since it is the pointwise maximum of a family of linear (thus convex) functions of $A$, it is convex.

8. **[Duality] Describe the dual of an SDP.**
   ??? success "Solution"
       $\max \sum b_i y_i$ subject to $\sum y_i A_i \preceq C$. The dual variables $y_i$ represent the shadow prices of the constraints.

9. **[Iterative] Mention a standard algorithm for solving large-scale SDPs.**
   ??? success "Solution"
       Primal-Dual Interior Point methods. They solve the KKT conditions iteratively using Newton's method on the central path.

10. **[Regularization] How does adding $\alpha \|X\|_F^2$ to an optimization problem affect the result?**
    ??? success "Solution"
        It acts as Tikhonov regularization, making the objective strictly convex and the solution unique. In matrix terms, it "pushes" the eigenvalues away from zero.

## Chapter Summary

This chapter explores the optimization landscape of linear operators:

1. **Variational Calculus**: Established eigenvalues as the extrema of quadratic forms, providing closed-form solutions to geometric problems.
2. **Spectral Interlacing**: Utilized the Courant-Fischer theorem to characterize the entire spectrum as a sequence of optimal subspaces.
3. **Trace Minimization**: Developed the calculus of trace functions for aligning and projecting high-dimensional data.
4. **Semidefinite Foundations**: Integrated matrix algebra with convex optimization to solve systems with positive definiteness constraints.
