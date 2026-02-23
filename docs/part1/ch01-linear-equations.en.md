# Chapter 01: Systems of Linear Equations

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Fields (Ch00) · Logic

**Chapter Outline**: Definition of Linear Systems → Augmented Matrix → Elementary Row Operations → Gaussian Elimination → Row Echelon Form (REF) and Reduced REF (RREF) → Consistency and Solvability → Free Variables and Rank → Homogeneous Systems → Inversion via Gauss-Jordan

**Extension**: Gaussian elimination is the fundamental algorithm of numerical linear algebra; its complexity $O(n^3)$ serves as the benchmark for matrix computation.

</div>

Systems of linear equations are the starting point of linear algebra. Every problem in the field—from finding eigenvalues to solving differential equations—eventually reduces to solving a system of the form $Ax = b$. The algorithm of **Gaussian elimination** provides a systematic way to solve these systems by transforming the augmented matrix into a simpler, row-echelon form.

---

## 01.1 Row Operations and Echelon Forms

!!! definition "Definition 01.1 (Elementary Row Operations)"
    1. **Swap**: Interchanging two rows.
    2. **Scale**: Multiplying a row by a non-zero scalar.
    3. **Pivot**: Adding a multiple of one row to another.
    These operations preserve the solution set of the system.

!!! theorem "Theorem 01.1 (Consistency)"
    A system $Ax = b$ is consistent (has at least one solution) if and only if the rank of the augmented matrix $[A | b]$ equals the rank of the coefficient matrix $A$.

---

## Exercises

1. **[Fundamentals] Solve the system $x+y=3, x-y=1$ using Gaussian elimination.**
   ??? success "Solution"
       Augmented matrix: $\begin{pmatrix} 1 & 1 & 3 \\ 1 & -1 & 1 \end{pmatrix} \to \begin{pmatrix} 1 & 1 & 3 \\ 0 & -2 & -2 \end{pmatrix} \to \begin{pmatrix} 1 & 1 & 3 \\ 0 & 1 & 1 \end{pmatrix} \to \begin{pmatrix} 1 & 0 & 2 \\ 0 & 1 & 1 \end{pmatrix}$. Solution is $x=2, y=1$.

2. **[RREF] What is the unique property of the Reduced Row Echelon Form (RREF)?**
   ??? success "Solution"
       While the steps of elimination may vary, the RREF of any matrix is unique. It reveals the pivot positions, the rank, and the basis for the column space directly.

3. **[Rank] If a $3 \times 5$ matrix has 3 pivots, what is its rank? How many free variables are in the system $Ax = 0$?**
   ??? success "Solution"
       The rank is 3 (number of pivots). The number of free variables is $n - \operatorname{rank}(A) = 5 - 3 = 2$.

4. **[Invertibility] How do you use Gauss-Jordan elimination to find $A^{-1}$?**
   ??? success "Solution"
       Perform row operations on the augmented matrix $[A | I]$ until the left side becomes $I$. The right side will then be $A^{-1}$. This succeeds iff $\operatorname{rank}(A) = n$.

5. **[Homogeneous] Show that a homogeneous system $Ax = 0$ with more variables than equations always has a non-trivial solution.**
   ??? success "Solution"
       If $n > m$, then $\operatorname{rank}(A) \le m < n$. The number of free variables $n - \operatorname{rank}(A) \ge 1$. Thus, there is at least one non-zero vector in the null space.

6. **[LU Link] Relate the steps of Gaussian elimination to the LU factorization.**
   ??? success "Solution"
       The row operations (without row swaps) can be recorded in a lower triangular matrix $L$. The final row echelon form is the upper triangular matrix $U$. Thus $A = LU$.

7. **[Consistency] For what value of $k$ is the system $\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix} x = \begin{pmatrix} 3 \\ k \end{pmatrix}$ consistent?**
   ??? success "Solution"
       Row 2 is twice Row 1. For consistency, $k$ must be twice 3, so $k=6$. Otherwise, the system represents parallel lines that never intersect.

8. **[Pivot] Why is choosing the largest entry as a pivot (Partial Pivoting) important in numerical computing?**
   ??? success "Solution"
       To reduce rounding errors. Dividing by a very small pivot can cause the coefficients of other rows to explode, leading to loss of precision in floating-point arithmetic.

9. **[Geometry] Describe the solution set of a $3 \times 3$ system with rank 2.**
   ??? success "Solution"
       The solution set is a line in $\mathbb{R}^3$ (assuming consistency). It represents the intersection of three planes that share a common line.

10. **[Complexity] Calculate the number of operations for Gaussian elimination on an $n \times n$ matrix.**
    ??? success "Solution"
        Approximately $2n^3/3$ floating-point operations (FLOPs). This cubic growth means that doubling the matrix size increases the computation time by eightfold.

## Chapter Summary

This chapter establishes the procedural engine of linear algebra:

1. **Algorithmic Consistency**: Defined elementary operations as the transformations that leave solution sets invariant.
2. **Structural Reduction**: Developed the row echelon form as the standard destination for system analysis.
3. **Existence and Uniqueness**: Formulated solvability criteria based on the rank of augmented matrices.
4. **Computational Baseline**: Linked the elimination process to matrix inversion and complexity analysis.
