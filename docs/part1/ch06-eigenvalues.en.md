# Chapter 06: Eigenvalues and Eigenvectors

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · Determinants (Ch3) · Linear Transformations (Ch5)

**Chapter Outline**: Eigenvalue Equation $Av = \lambda v$ → Characteristic Polynomial → Algebraic and Geometric Multiplicity → Conditions for Matrix Diagonalization → Similarity Transformation $P^{-1}AP = D$ → Spectral Theorem for Symmetric Matrices → Cayley-Hamilton Theorem

**Extension**: Eigenvalues reveal the pure scaling behavior of a linear operator in certain specific directions and are the golden key to analyzing the stability of dynamical systems.

</div>

Eigenvalue theory is one of the most profound parts of linear algebra. It reduces complex matrix operations to scalar multiplications. By finding the "natural frequencies" (eigenvalues) and "mode shapes" (eigenvectors) of a matrix, we can deconstruct the core structure of an operator.

---

## 06.1 Core Definitions

!!! definition "Definition 06.1 (Eigenvalues and Eigenvectors)"
    For an $n \times n$ matrix $A$, if there exists a non-zero vector $v$ and a scalar $\lambda$ such that:
    $$Av = \lambda v$$
    then $\lambda$ is an **eigenvalue** of $A$ and $v$ is an **eigenvector** corresponding to $\lambda$.

!!! theorem "Theorem 06.3 (Relation to Trace and Determinant)"
    1. $\sum \lambda_i = \operatorname{tr}(A)$
    2. $\prod \lambda_i = \det(A)$

---

## Exercises

1. **[Basic Calculation] Calculate the eigenvalues of $A = \begin{pmatrix} 4 & 1 \\ 2 & 3 \end{pmatrix}$.**
   ??? success "Solution"
       Characteristic equation: $\det(\lambda I - A) = (\lambda-4)(\lambda-3) - 2 = \lambda^2 - 7\lambda + 10 = 0$.
       Solving gives $\lambda_1 = 5, \lambda_2 = 2$.

2. **[Eigenvectors] Find the eigenvector corresponding to $\lambda=5$ for $A = \begin{pmatrix} 4 & 1 \\ 2 & 3 \end{pmatrix}$.**
   ??? success "Solution"
       Solve $(5I-A)v = 0 \implies \begin{pmatrix} 1 & -1 \\ -2 & 2 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = 0$.
       Thus $x_1 = x_2$. The eigenvector is $c \begin{pmatrix} 1 \\ 1 \end{pmatrix}$ (for $c \neq 0$).

3. **[Diagonalization] Determine if $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$ is diagonalizable.**
   ??? success "Solution"
       Eigenvalues are $\lambda = 1, 1$ (algebraic multiplicity 2).
       Solving for eigenvectors: $(I-A)v = 0 \implies \begin{pmatrix} 0 & -1 \\ 0 & 0 \end{pmatrix} v = 0 \implies x_2 = 0$.
       The dimension of the eigenspace (geometric multiplicity) is 1. Since geometric multiplicity < algebraic multiplicity, it is not diagonalizable.

4. **[Trace and Det] If the eigenvalues of a $3 \times 3$ matrix are 1, 2, 3, find $\operatorname{tr}(A)$ and $\det(A)$.**
   ??? success "Solution"
       $\operatorname{tr}(A) = 1+2+3 = 6$.
       $\det(A) = 1 \cdot 2 \cdot 3 = 6$.

5. **[Spectral Theorem Intro] Prove: The eigenvalues of a real symmetric matrix are always real.**
   ??? success "Solution"
       Let $Av = \lambda v$. Then $\bar{v}^T Av = \lambda \bar{v}^T v$.
       Taking the conjugate transpose and using $A^T=A$ gives $\bar{v}^T A v = \bar{\lambda} \bar{v}^T v$.
       Thus $\lambda \|v\|^2 = \bar{\lambda} \|v\|^2$. Since $v \neq 0$, we must have $\lambda = \bar{\lambda}$, so $\lambda$ is real.

6. **[Powers] If $\lambda$ is an eigenvalue of $A$, prove $\lambda^k$ is an eigenvalue of $A^k$.**
   ??? success "Solution"
       $Av = \lambda v \implies A(Av) = A(\lambda v) = \lambda (Av) = \lambda^2 v$.
       By induction, $A^k v = \lambda^k v$.

7. **[Cayley-Hamilton] Verify that $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$ satisfies its characteristic equation.**
   ??? success "Solution"
       Characteristic polynomial $p(\lambda) = (\lambda-1)(\lambda-2) = \lambda^2 - 3\lambda + 2$.
       $p(A) = A^2 - 3A + 2I = \begin{pmatrix} 1 & 0 \\ 0 & 4 \end{pmatrix} - \begin{pmatrix} 3 & 0 \\ 0 & 6 \end{pmatrix} + \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$.

8. **[Inverses] If $A$ is invertible with eigenvalue $\lambda$, prove $1/\lambda$ is an eigenvalue of $A^{-1}$.**
   ??? success "Solution"
       $Av = \lambda v \implies v = A^{-1}(\lambda v) = \lambda A^{-1} v$.
       Since $A$ is invertible, $\lambda \neq 0$. Dividing by $\lambda$ gives $A^{-1} v = (1/\lambda) v$.

9. **[Similarity] Prove that similar matrices have the same eigenvalues.**
   ??? success "Solution"
       $\det(\lambda I - P^{-1}AP) = \det(P^{-1}(\lambda I - A)P) = \det(P^{-1})\det(\lambda I - A)\det(P) = \det(\lambda I - A)$. The characteristic polynomials are identical, hence eigenvalues are the same.

10. **[Algebraic Multiplicity] What is the algebraic multiplicity of eigenvalue 2 for $A = \begin{pmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 2 \end{pmatrix}$?**
    ??? success "Solution"
        The characteristic polynomial is $(\lambda-2)^3$. Thus, the algebraic multiplicity is 3.

## Chapter Summary

Eigenvalue theory is the soul of linear algebra:

1. **Decoupling**: Diagonalization is essentially breaking a complex coupled system into independent linear components.
2. **Multiplicity Balance**: The matching of algebraic and geometric multiplicities is the key to matrix simplification.
3. **Structural Conservation**: Trace and determinant, as functions of eigenvalues, are the cornerstones of invariant theory.
