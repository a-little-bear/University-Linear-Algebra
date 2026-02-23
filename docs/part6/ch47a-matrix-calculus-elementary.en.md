# Chapter 47A: Matrix Calculus - Elementary Theory

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Multivariable Calculus · Kronecker Products (Ch19)

**Chapter Outline**: Gradient of Scalars with Respect to Matrices → Derivative of Matrices with Respect to Scalars → Core Derivative Formulas (Trace, Determinant, Inverse) → Matrix Form of the Chain Rule → The Commutation Matrix $K_{mn}$ → The Duplication Matrix $D_n$ → Differentiation under Symmetry Constraints → Applications: Least Squares Optimization and Maximum Likelihood Estimation

**Extension**: Matrix calculus is the mathematical foundation for backpropagation in modern deep learning; it integrates scattered multivariable partial derivatives into compact operator forms, making the optimization of complex models (like covariance structures) possible.

</div>

Matrix calculus extends the concept of derivatives from classical calculus to matrix variables. When dealing with optimization problems involving thousands of variables (such as neural network training), explicitly writing out every partial derivative is impractical. Through matrix calculus, we arrange partial derivatives into gradient matrices isomorphic to the original matrices, utilizing concise algebraic notation to complete the differentiation of all variables at once.

---

## 47A.1 Gradients and Jacobian Matrices

!!! definition "Definition 47A.1 (Gradient of a Scalar w.r.t. a Matrix)"
    Let $f(X)$ be a scalar function of an $m 	imes n$ matrix $X$. Its gradient is defined as a matrix of the same dimensions:
    $$
abla_X f(X) = \frac{\partial f}{\partial X} = \left[ \frac{\partial f}{\partial x_{ij}} ight]$$

!!! definition "Definition 47A.2 (Derivative of a Matrix w.r.t. a Scalar)"
    Let $X(t)$ be a matrix-valued function of a scalar $t$. Its derivative is defined as:
    $$\frac{dX}{dt} = \left[ \frac{dx_{ij}}{dt} ight]$$

---

## 47A.2 Core Differentiation Formulas

!!! theorem "Theorem 47A.1 (Derivatives of Trace and Determinant)"
    1.  **Linearity of Trace**: $
abla_X \operatorname{tr}(AX) = A^T$
    2.  **Quadratic Form**: $
abla_x (\mathbf{x}^T A \mathbf{x}) = (A + A^T) \mathbf{x}$
    3.  **Determinant (Jacobi’s Formula)**: $
abla_X \det(X) = \det(X) X^{-T}$
    4.  **Log-Determinant**: $
abla_X \ln \det(X) = X^{-T}$
    5.  **Inverse Matrix**: $\frac{d}{dt} (X^{-1}) = -X^{-1} \frac{dX}{dt} X^{-1}$

---

## 47A.3 Structural Matrices

!!! technique "Commutation Matrix $K_{mn}$"
    The **Commutation Matrix** is the unique permutation matrix satisfying $\operatorname{vec}(A^T) = K_{mn} \operatorname{vec}(A)$. It acts as an "index swapper" in the chain rule involving transpose terms.

!!! technique "Duplication Matrix $D_n$"
    When differentiating with respect to **Symmetric Matrices**, the variables $x_{ij} = x_{ji}$ are not independent. The duplication matrix maps unique lower-triangular elements to the full vectorized symmetric matrix, correcting the gradient.

---

## Exercises

1.  **[Calculation] Find $
abla_x (\|\mathbf{x}\|_2^2)$.**
??? success "Solution"
     $
abla_x (\mathbf{x}^T \mathbf{x}) = 2\mathbf{x}$.

2.  **[Trace] Find $
abla_X \operatorname{tr}(X^T A X)$.**
??? success "Solution"
     $(A + A^T) X$.


****
??? success "Solution"
     $\operatorname{tr}(-X^{-1} \dot{X} X^{-1}) = -\operatorname{tr}(X^{-2} \dot{X})$.

4.  **[Determinant] If $A$ is a constant matrix, find $
abla_X \det(AX)$.**
??? success "Solution"
     $\det(A) 
abla_X \det(X) = \det(A)\det(X) X^{-T} = \det(AX) X^{-T}$.


****
??? success "Solution"
     This follows directly from the vectorization identity of the Kronecker product.


****
??? success "Solution"
     Because $x_{ij}$ and $x_{ji}$ are the same variable, the result is $A + A^T - \operatorname{diag}(A)$.


****
??? success "Solution"
     $
abla_x y = \left( \frac{\partial u}{\partial x} ight)^T 
abla_u f$.

8.  **[Frobenius] Find $
abla_X (\|X\|_F^2)$.**
??? success "Solution"
     $
abla_X \operatorname{tr}(X^T X) = 2X$.


****
??? success "Solution"
     It is simply the matrix $A$.

****
??? success "Solution"
    
abla_x (y-Ax)^T(y-Ax) = -2A^T(y-Ax) = 0 \implies A^T A x = A^T y$ (the Normal Equations).

## Chapter Summary

Matrix calculus enables the leap from scalar partial derivatives to operator gradients:


****: Through gradient matrices and Jacobians, complex multivariate changes are condensed into single algebraic terms, greatly simplifying the derivation of high-dimensional optimization problems.

****: Core derivative formulas (such as those for the determinant and inverse) reveal the non-linear sensitivity of global matrix properties to local entries, serving as the bedrock of stability analysis.

****: The introduction of commutation and duplication matrices provides standard algebraic compensation for inherent matrix symmetry and index permutations, ensuring the rigor of the differentiation process.
