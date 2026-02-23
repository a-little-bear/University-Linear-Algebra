# Chapter 06: Eigenvalues and Eigenvectors

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Determinants (Ch3) · Linear Transformations (Ch5) · Polynomials (Ch00)

**Chapter Outline**: The Eigenvalue Equation $Av = \lambda v$ → Characteristic Polynomial $\det(\lambda I - A) = 0$ → Algebraic and Geometric Multiplicities → Diagonalization ($A = PDP^{-1}$) → Powers of Matrices → Eigenspaces → Defective Matrices → Trace and Determinant via Eigenvalues → Cayley-Hamilton Theorem → Gershgorin Disc Theorem

**Extension**: Eigenvalues reveal the "natural frequencies" and scaling factors of an operator; diagonalization is the ultimate tool for solving linear difference and differential equations.

</div>

Eigenvalues and eigenvectors capture the invariant directions of a linear transformation. When a matrix $A$ acts on an eigenvector $v$, it simply scales it by a factor $\lambda$ without changing its direction. This "spectral" decomposition allows us to decouple complex systems of equations and simplify the calculation of matrix powers and exponentials. This chapter details the process of finding eigenvalues via the characteristic polynomial and establishes the criteria for matrix **diagonalization**.

---

## 06.1 Eigenvalue Theory and Diagonalization

!!! definition "Definition 06.1 (Eigenvalue Equation)"
    A scalar $\lambda$ is an **eigenvalue** of $A$ and $v \neq 0$ is the corresponding **eigenvector** if:
    $$Av = \lambda v$$
    This is equivalent to $(\lambda I - A)v = 0$, meaning $(\lambda I - A)$ must be singular.

!!! theorem "Theorem 06.1 (Diagonalization Criterion)"
    An $n \times n$ matrix $A$ is diagonalizable if and only if it has $n$ linearly independent eigenvectors. This occurs iff the geometric multiplicity equals the algebraic multiplicity for every eigenvalue.

---

## Exercises

1. **[Fundamentals] Find the eigenvalues of $A = \begin{pmatrix} 4 & 1 \\ 2 & 3 \end{pmatrix}$.**
   ??? success "Solution"
       $\det(\lambda I - A) = (\lambda-4)(\lambda-3) - 2 = \lambda^2 - 7\lambda + 10 = (\lambda-2)(\lambda-5) = 0$. Eigenvalues are $\lambda_1 = 2, \lambda_2 = 5$.

2. **[Eigenvectors] Find an eigenvector for $\lambda = 2$ from the previous matrix.**
   ??? success "Solution"
       $(2I - A)v = \begin{pmatrix} -2 & -1 \\ -2 & -1 \end{pmatrix} \begin{pmatrix} x \\ y \end{pmatrix} = 0$. This implies $2x+y=0$. An eigenvector is $v = (1, -2)^T$.

3. **[Diagonalization] Diagonalize $A = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}$.**
   ??? success "Solution"
       Eigenvalues are 1 and 2. Eigenvectors are $v_1 = (1, 0)^T$ and $v_2 = (1, 1)^T$. Let $P = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$ and $D = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$. Then $A = PDP^{-1}$.

4. **[Powers] Calculate $A^{10}$ for $A = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}$.**
   ??? success "Solution"
       $A^{10} = P D^{10} P^{-1} = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} 1 & 0 \\ 0 & 2^{10} \end{pmatrix} \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 2^{10}-1 \\ 0 & 2^{10} \end{pmatrix}$.

5. **[Trace and Det] Relate $\operatorname{tr}(A)$ and $\det A$ to eigenvalues.**
   ??? success "Solution"
       $\operatorname{tr}(A) = \sum \lambda_i$ and $\det A = \prod \lambda_i$. This provides a quick check for eigenvalue calculations.

6. **[Cayley-Hamilton] What does the Cayley-Hamilton theorem state?**
   ??? success "Solution"
       Every square matrix satisfies its own characteristic polynomial: $p(A) = 0$. This allows for expressing $A^n$ and $A^{-1}$ as polynomials in $A$ of degree $< n$.

7. **[Defect] Can $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ be diagonalized?**
   ??? success "Solution"
       No. Eigenvalues are 0 (algebraic multiplicity 2). $(0I-A) = \begin{pmatrix} 0 & -1 \\ 0 & 0 \end{pmatrix}$ has a null space of dimension 1. Geometric multiplicity (1) < Algebraic multiplicity (2). The matrix is **defective**.

8. **[Invertibility] If $\lambda=0$ is an eigenvalue of $A$, is $A$ invertible?**
   ??? success "Solution"
       No. $\det A = \prod \lambda_i$. If one eigenvalue is 0, then $\det A = 0$, so $A$ is singular.

9. **[Symmetry] Why are real symmetric matrices always diagonalizable?**
   ??? success "Solution"
       The Spectral Theorem ensures that for symmetric $A$, eigenvectors corresponding to distinct eigenvalues are orthogonal, and there is always a full set of $n$ orthonormal eigenvectors.

10. **[Gershgorin] Use the Gershgorin Disc Theorem to bound the eigenvalues of $\begin{pmatrix} 10 & 1 \\ 2 & 5 \end{pmatrix}$.**
    ??? success "Solution"
        Discs are $D_1: |z-10| \le 1$ and $D_2: |z-5| \le 2$. All eigenvalues lie in the union of these two circles in the complex plane.

## Chapter Summary

This chapter establishes the spectral analysis of linear operators:

1. **Invariant Subspaces**: Defined eigenvalues and eigenvectors as the fundamental descriptors of scaling and direction.
2. **Structural Simplification**: Developed diagonalization as the process of converting a complex operator into a simple coordinate-wise scaling.
3. **Multiplicity Logic**: Distinguished between algebraic and geometric aspects of eigenvalues to identify defective matrices.
4. **Polynomial Constraint**: Leveraged the Cayley-Hamilton theorem to link matrix powers to the characteristic polynomial.
