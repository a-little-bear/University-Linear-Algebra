# Chapter 39: Totally Positive Matrices

<div class="context-flow" markdown>

**Prerequisites**: Determinants (Ch3) · Matrix Factorization (Ch10) · Eigenvalues (Ch6) · Sign Patterns (Ch65A)

**Chapter Outline**: Definition of Total Positivity (TP) and Total Nonnegativity (TN) → Criteria for TP → Properties of Minors → Variation Diminishing Property → Eigenvalue Properties (Real and Distinct) → Oscillation Matrices → LU Factorization (Bidiagonal Factorization) → Application in Spline Theory and Statistics

**Extension**: Totally positive matrices are the foundation of spline interpolation and the analysis of variation-diminishing linear maps.

</div>

A matrix is **totally positive** (TP) if *every* one of its minors (not just principal ones) is strictly positive. This is an extremely restrictive condition that leads to remarkable properties: the eigenvalues are real, positive, and distinct, and the associated linear transformation is **variation-diminishing** (it never increases the number of sign changes in a vector). TP matrices serve as the discrete analog of kernels for many integral equations in analysis.

---

## 39.1 Definitions and Structural Properties

!!! definition "Definition 39.1 (Totally Positive and Nonnegative)"
    A matrix $A$ is **totally positive** (TP) if all its minors are strictly positive.
    A matrix $A$ is **totally nonnegative** (TN) if all its minors are non-negative.

!!! theorem "Theorem 39.1 (Eigenvalue Rigidity)"
    If $A$ is an $n 	imes n$ totally positive matrix, then its eigenvalues are real, positive, and simple (distinct):
    $$\lambda_1 > \lambda_2 > \dots > \lambda_n > 0$$

---

## Exercises

1. **[Fundamentals] Is the Vandermonde matrix $V(x_1, \dots, x_n)$ totally positive for $0 < x_1 < \dots < x_n$?**
   ??? success "Solution"
       Yes. Every minor of a Vandermonde matrix with positive increasing nodes is a Schur polynomial times a product of differences, which is strictly positive.

2. **[Criteria] How many minors do you need to check to verify if an $n 	imes n$ matrix is TP?**
   ??? success "Solution"
       While there are $\binom{2n}{n}-1$ total minors, checking only the $n^2$ **contiguous** minors (minors formed by adjacent rows and columns) is sufficient to prove total positivity.

3. **[Variation Diminishing] Define the variation diminishing property of a TP matrix.**
   ??? success "Solution"
       Let $v(x)$ be the number of sign changes in vector $x$. If $A$ is TP, then $v(Ax) \le v(x)$. The transformation $A$ "smoothes" the signal, never creating new oscillations.

4. **[Cauchy-Binet] Prove that the set of TN matrices is closed under multiplication.**
   ??? success "Solution"
       By the Cauchy-Binet formula, a minor of $AB$ is a sum of products of minors of $A$ and $B$. Since all minors of $A$ and $B$ are non-negative, their sum is non-negative. Thus $AB$ is TN.

5. **[Oscillation Matrices] What is an oscillation matrix?**
   ??? success "Solution"
       A TN matrix $A$ is an oscillation matrix if there exists some power $k$ such that $A^k$ is TP. This class shares the distinct positive eigenvalue property of TP matrices.

6. **[Factorization] Describe the bidiagonal factorization of a TN matrix.**
   ??? success "Solution"
       Every TN matrix $A$ can be factored into a product of non-negative elementary bidiagonal matrices. This provides a way to parameterize TN matrices using only non-negative numbers.

7. **[Pascal Matrix] Show that the Pascal matrix $P_{ij} = \binom{i+j}{i}$ is totally positive.**
   ??? success "Solution"
       The Pascal matrix can be factored as $L \cdot U$, where $L$ and $U$ are lower and upper triangular Pascal matrices. Each factor is TN, and their product is TP.

8. **[Inversion] Is the inverse of a TP matrix also TP?**
   ??? success "Solution"
       No. The inverse of a TP matrix has a checkerboard sign pattern ($(-1)^{i+j} (A^{-1})_{ij} > 0$). Such matrices are called **strictly sign-regular**.

9. **[Splines] Why are B-splines related to total positivity?**
   ??? success "Solution"
       The collocation matrix of B-spline basis functions is totally positive. This ensures that spline approximation is stable and that the spline curve follows the "shape" of the control points (variation diminishing).

10. **[Compound Matrices] Relate TP matrices to the eigenvalues of compound matrices.**
    ??? success "Solution"
        A matrix $A$ is TP iff all its $k$-th order compound matrices $C_k(A)$ (matrices of $k 	imes k$ minors) are positive matrices. This allows the use of Perron-Frobenius theory on the compound matrices to prove that the eigenvalues of $A$ are distinct.

## Chapter Summary

This chapter explores the rigid structure of matrices with all positive minors:

1. **Global Positivity**: Defined TP and TN matrices as operators where every sub-volume is positively oriented.
2. **Spectral Separation**: Established that total positivity forces eigenvalues to be real, positive, and strictly separated.
3. **Shape Preservation**: Analyzed the variation-diminishing property, linking matrix algebra to the geometry of curves and splines.
4. **Factorization Logic**: Demonstrated that TP matrices are generated by the product of non-negative bidiagonal transformations.
