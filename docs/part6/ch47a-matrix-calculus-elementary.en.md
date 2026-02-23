# Chapter 47A: Matrix Calculus: Element-wise and Structural

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Multivariable Calculus · Kronecker Product (Ch19) · Matrix Norms (Ch15)

**Chapter Outline**: Differential of a Matrix Function → Gradient, Jacobian, and Hessian → Trace Trick for Differentials → Derivatives of Inverses and Determinants → Chain Rule for Matrices → Vec Operator and Kronecker Products → Optimization on Matrix Manifolds → Application in Machine Learning (Backpropagation)

**Extension**: Matrix calculus is the mathematical engine behind modern deep learning (optimizing weights in neural networks) and Gaussian process regression.

</div>

Matrix calculus extends the rules of differentiation to functions of matrices. Instead of dealing with $n^2$ individual partial derivatives, we treat the matrix as a single variable. The core technique is the **differential method**: using the identity $df = \operatorname{tr}(G^T dX)$ to identify the gradient $G$. This approach is "layout-agnostic" and simplifies complex expressions involving inverses, determinants, and traces.

---

## 47A.1 Differentials and Gradients

!!! definition "Definition 47A.1 (Gradient of a Scalar Function)"
    For a scalar function $f(X)$ where $X \in M_{m \times n}$, the **gradient** $\nabla_X f$ is the matrix of the same dimension:
    $$(\nabla_X f)_{ij} = \frac{\partial f}{\partial X_{ij}}$$
    The differential is related to the gradient by: $df = \langle \nabla_X f, dX \rangle = \operatorname{tr}((\nabla_X f)^T dX)$.

!!! theorem "Theorem 47A.1 (The "Big Three" Differentials)"
    1. $d(AXB) = A(dX)B$
    2. $d(X^{-1}) = -X^{-1}(dX)X^{-1}$
    3. $d(\det X) = \det X \cdot \operatorname{tr}(X^{-1} dX)$

---

## Exercises

1. **[Fundamentals] Find the gradient of $f(x) = x^T A x$.**
   ??? success "Solution"
       $df = d(x^T A x) = (dx)^T A x + x^T A (dx) = x^T A^T dx + x^T A dx = x^T (A + A^T) dx$. Thus $\nabla_x f = (A + A^T) x$. If $A$ is symmetric, $\nabla_x f = 2Ax$.

2. **[Trace Trick] Find the gradient of $f(X) = \operatorname{tr}(AX)$.**
   ??? success "Solution"
       $df = \operatorname{tr}(A dX)$. By the definition $df = \operatorname{tr}(G^T dX)$, we have $G^T = A$, so $\nabla_X f = A^T$.

3. **[Inverse] Compute the derivative of $f(X) = \operatorname{tr}(X^{-1}A)$.**
   ??? success "Solution"
       $df = \operatorname{tr}(d(X^{-1})A) = \operatorname{tr}(-X^{-1}(dX)X^{-1}A) = \operatorname{tr}(-X^{-1}AX^{-1} dX)$. Thus $\nabla_X f = -(X^{-1}AX^{-1})^T = -X^{-T}A^TX^{-T}$.

4. **[Log-Det] Derive the gradient of $f(X) = \log\det X$ for $X \succ 0$.**
   ??? success "Solution"
       $df = d(\log\det X) = \frac{1}{\det X} d(\det X) = \frac{1}{\det X} (\det X \operatorname{tr}(X^{-1} dX)) = \operatorname{tr}(X^{-1} dX)$. Thus $\nabla_X f = X^{-T} = X^{-1}$ (since $X$ is symmetric).

5. **[Chain Rule] Find the gradient of $f(X) = \|AX-B\|_F^2$.**
   ??? success "Solution"
       Let $R = AX-B$. $f = \operatorname{tr}(R^T R)$. $df = 2 \operatorname{tr}(R^T dR) = 2 \operatorname{tr}(R^T A dX) = 2 \operatorname{tr}(R^T A dX)$. Thus $\nabla_X f = (2 R^T A)^T = 2 A^T (AX-B)$. Setting this to zero gives the normal equations for least squares.

6. **[Vec Operator] Use the identity $\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$ to find the Jacobian of $f(X) = AXB$.**
   ??? success "Solution"
       $\frac{\partial \operatorname{vec}(AXB)}{\partial \operatorname{vec}(X)} = B^T \otimes A$. This allows for converting matrix calculus into vector calculus operations.

7. **[Product Rule] Show that $d(XY) = (dX)Y + X(dY)$.**
   ??? success "Solution"
       By the incremental definition: $(X+dX)(Y+dY) - XY = XY + (dX)Y + X(dY) + (dX)(dY) - XY$. Dropping the second-order term $(dX)(dY)$ gives the result.

8. **[Hessian] Define the Hessian of a scalar function $f(X)$.**
   ??? success "Solution"
       The Hessian is the second derivative matrix. In matrix calculus, it is often represented as a bilinear operator such that $d^2 f = \operatorname{vec}(dX)^T \mathbf{H} \operatorname{vec}(dX)$.

9. **[Eigenvalues] What is the differential of the $i$-th eigenvalue $d\lambda_i$ for a symmetric matrix?**
   ??? success "Solution"
       $d\lambda_i = v_i^T (dA) v_i$, where $v_i$ is the corresponding unit eigenvector. This result is essential for sensitivity analysis in vibration theory and statistics.

10. **[Complexity] Why is the differential method superior to the element-by-element method?**
    ??? success "Solution"
        It treats the matrix as a single algebraic object, preserving structural identities (like the cyclic trace property). This leads to coordinate-free results that are easier to verify and implement in high-level programming languages.

## Chapter Summary

This chapter establishes the calculus of matrix-valued variables:

1. **Differential Calculus**: Developed the differential method as the standard tool for identifying gradients and Jacobians.
2. **Structural Identities**: Derived the derivatives of fundamental matrix operators, including inverses and determinants.
3. **Kronecker Bridge**: Utilized the vec-Kronecker identity to link matrix calculus to vector-based optimization.
4. **Optimization Engine**: Demonstrated the application of matrix derivatives in least-squares, likelihood maximization, and neural network training.
