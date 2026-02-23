# Chapter 25: Applications of Linear Algebra in Optimization

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Matrix Calculus (Ch47) · SVD (Ch11) · Matrix Manifolds (Ch24)

**Chapter Outline**: Modeling Optimization via Linear Algebra → Unconstrained Optimization (Gradient Descent, Newton's Method) → Linear Programming (LP) and Simplex Method → Quadratic Programming (QP) → Semidefinite Programming (SDP) → Interior Point Methods → Lagrange Multipliers and Duality → Compressed Sensing and $L_1$ Minimization → ADMM

**Extension**: Optimization is the largest "consumer" of linear algebra theory; from determining search directions to proving convergence, linear algebra provides all the analytical tools.

</div>

Optimization is the process of finding extremum points of functions under constraints. In this field, matrices are not just containers for coefficients, but core operators that define geometric curvature, gradient directions, and variable coupling. Almost all modern complex optimization algorithms boil down to repeatedly solving linear systems or eigenvalue problems.

---

## 25.1 Core Models and Duality Theory

!!! definition "Definition 25.1 (Semidefinite Programming - SDP)"
    SDP involves minimizing a linear function subject to linear constraints and the requirement that the variable matrix $X$ is positive semi-definite:
    $$\min \operatorname{tr}(CX), 	ext{ s.t. } \operatorname{tr}(A_i X) = b_i, X \succeq 0$$

!!! theorem "Theorem 25.3 (Strong Duality)"
    If both the primal and dual problems are strictly feasible, the optimal value of the primal problem equals the optimal value of the dual. This provides the theoretical basis for verifying optimality.

---

## Exercises

1. **[Newton's Method] Why does Newton's method converge faster than gradient descent? Explain using linear algebra.**
   ??? success "Solution"
       Gradient descent uses only first-order information (direction), while Newton's method utilizes second-order information (the inverse Hessian $H^{-1}$). $H^{-1}$ acts as a "preconditioner" that corrects the gradient direction based on local curvature, allowing the search to point directly toward the minimum on ill-conditioned (long and thin) surfaces.

2. **[Quadratic Programming] Write the formula for the minimum point of an unconstrained quadratic form $f(x) = \frac{1}{2}x^T Ax - b^T x$.**
   ??? success "Solution"
       Differentiating gives $
abla f = Ax - b = 0$. If $A \succ 0$, the minimum is $x^* = A^{-1} b$. This shows that solving a linear system is equivalent to minimizing a quadratic energy function.

3. **[Lagrange Multipliers] Solve $\min \|x\|^2$ subject to $Ax=b$.**
   ??? success "Solution"
       Construct $L(x, \lambda) = x^T x + \lambda^T(Ax-b)$.
       $
abla_x L = 2x + A^T \lambda = 0 \implies x = -0.5 A^T \lambda$.
       Substitute into constraint: $A(-0.5 A^T \lambda) = b \implies \lambda = -2(AA^T)^{-1}b$.
       Thus $x^* = A^T(AA^T)^{-1}b$. This is the **minimum norm solution** (using the generalized inverse).

4. **[SDP Application] Why is SDP useful in combinatorial optimization?**
   ??? success "Solution"
       Many hard discrete constraints (like $x_i \in \{-1, 1\}$) can be relaxed into matrix rank constraints, which are then further relaxed into convex semi-definite constraints $X \succeq 0$. SDP provides powerful polynomial-time lower bound estimates (e.g., for the Max-Cut problem).

5. **[KKT Conditions] Write the complementary slackness condition in the KKT conditions for constrained optimization.**
   ??? success "Solution"
       $\lambda_i g_i(x^*) = 0$. This means that either the constraint is inactive ($\lambda_i=0$) or the solution lies exactly on the boundary ($g_i(x^*)=0$).

6. **[Condition Number] Why is it difficult for optimization algorithms to converge when the matrix $A$ has a large condition number?**
   ??? success "Solution"
       A large condition number means the level sets of the objective function are extremely flat ellipsoids. The gradient direction is nearly orthogonal to the direction pointing toward the minimum, causing the algorithm to oscillate across the "valley" walls, making the step size difficult to determine.

7. **[ADMM] Briefly describe the core idea of ADMM.**
   ??? success "Solution"
       ADMM decomposes the augmented Lagrangian function to split a large coupled problem into multiple small, easy-to-solve sub-problems (alternatingly optimizing variables). It is particularly suited for large-scale, distributed optimization with non-smooth terms (like the $L_1$ norm).

8. **[Compressed Sensing] Why does $L_1$ minimization produce sparse solutions?**
   ??? success "Solution"
       The level sets of the $L_1$ norm (diamonds) have "corners" on the coordinate axes. When the optimization process looks for a tangent point with the linear constraint plane, it is highly likely to hit a corner on an axis, forcing many components to zero.

9. **[Diagonal Dominance App] Why is $A + \lambda I$ often used to replace an unstable $A$ in numerical optimization?**
   ??? success "Solution"
       This is called **Tikhonov Regularization** or ridge regression. Adding $\lambda I$ ($\lambda > 0$) forces the matrix to be strictly diagonally dominant, improving its condition number and ensuring the inverse exists and is numerically stable.

10. **[SVD and PCA] Prove that Principal Component Analysis (PCA) is essentially a trace maximization problem.**
    ??? success "Solution"
        PCA aims to find a projection matrix $P$ that maximizes $\operatorname{tr}(P^T \Sigma P)$. According to the Ky Fan Maximum Principle, the optimal $P$ is composed of the eigenvectors corresponding to the $k$ largest eigenvalues of the covariance matrix $\Sigma$.

## Chapter Summary

Optimization is the "action plan" of linear algebra:

1. **Duality of Equations and Extrema**: Solving equations is finding extrema; finding extrema is solving equations.
2. **The Language of Curvature**: The spectral properties of the Hessian matrix determine search efficiency and solution stability.
3. **The Art of Relaxation**: Embedding non-convex, discrete problems into convex semi-definite spaces through matrixization is the essence of modern optimization.
