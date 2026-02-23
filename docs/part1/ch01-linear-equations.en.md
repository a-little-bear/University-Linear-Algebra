# Chapter 01: Linear Equations

<div class="context-flow" markdown>

**Prerequisites**: Polynomial Algebra (Ch00)

**Chapter Outline**: Definition of Linear Systems → Matrix Representation (Augmented Matrices) → Elementary Row Operations → Row Echelon Form (REF) & Reduced Row Echelon Form (RREF) → Gaussian & Gauss-Jordan Elimination → Existence and Uniqueness Theorems (Rank-Presence) → Homogeneous & Non-homogeneous Systems → Geometric Interpretation (Intersections of Hyperplanes) → Introduction to Numerical Stability

**Extension**: Solving linear systems is the underlying engine of nearly all numerical computing (Finite Elements, Optimization, Machine Learning); the structure of its solution space leads directly to the definition of Vector Spaces (Ch04).

</div>

Systems of linear equations are the logical starting point of linear algebra. From ancient algorithmic texts to modern supercomputers, solving $Ax = b$ remains the most central task in scientific computing. This chapter establishes the standard framework for handling linear systems from two dimensions: algorithmic (elimination methods) and theoretical (the structure of solutions).

---

## 01.1 Linear Systems and Matrix Representation

!!! definition "Definition 01.1 (System of Linear Equations)"
    A system of $m$ linear equations in $n$ variables is typically written as:
    $$\begin{cases} a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n = b_1 \\ a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n = b_2 \\ \vdots \\ a_{m1}x_1 + a_{m2}x_2 + \cdots + a_{mn}x_n = b_m \end{cases}$$
    where $a_{ij}$ are coefficients, $b_i$ are constants, and $x_j$ are variables.

!!! definition "Definition 01.2 (Augmented Matrix)"
    By combining the coefficients and constants into a single matrix, we obtain the $m \times (n+1)$ **augmented matrix**:
    $$\tilde{A} = [A | \mathbf{b}] = \begin{pmatrix} a_{11} & a_{12} & \cdots & a_{1n} & | & b_1 \\ a_{21} & a_{22} & \cdots & a_{2n} & | & b_2 \\ \vdots & \vdots & \ddots & \vdots & | & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} & | & b_m \end{pmatrix}$$

---

## 01.2 Elementary Row Operations and Echelon Forms

!!! definition "Definition 01.3 (Elementary Row Operations)"
    1.  **Scaling**: Multiply a row by a non-zero constant $k$.
    2.  **Swapping**: Interchange two rows.
    3.  **Replacement**: Add a multiple of one row to another row.
    These operations **do not change** the solution set of the system.

!!! definition "Definition 01.4 (Reduced Row Echelon Form - RREF)"
    A matrix is in **Reduced Row Echelon Form** (RREF) if:
    1.  All non-zero rows are above any rows of all zeros.
    2.  The leading entry (pivot) of each non-zero row is 1.
    3.  Each leading 1 is the only non-zero entry in its column.
    4.  The leading 1 of a row is to the right of the leading 1 of the row above it.

---

## 01.3 Gaussian Elimination

!!! algorithm "Algorithm 01.1 (Gauss-Jordan Elimination Steps)"
    1.  **Forward Phase**: Use elementary operations to transform the augmented matrix into Row Echelon Form (REF).
    2.  **Backward Phase**: Scaling pivots to 1 and clearing entries above them to reach RREF.
    3.  **Read Solution**: Write the general solution directly from the RREF.

---

## 01.4 Existence and Uniqueness

!!! theorem "Theorem 01.1 (Existence and Uniqueness Theorem)"
    Let $A$ be an $m \times n$ matrix, and let $\operatorname{rank}(A)$ be its rank.
    1.  **No Solution (Inconsistent)**: $\operatorname{rank}(A) < \operatorname{rank}([A|\mathbf{b}])$. This happens if a row like $[0 \ 0 \ \cdots \ 0 \ | \ d]$ ($d \neq 0$) appears in RREF.
    2.  **Unique Solution**: $\operatorname{rank}(A) = \operatorname{rank}([A|\mathbf{b}]) = n$ (number of variables).
    3.  **Infinitely Many Solutions**: $\operatorname{rank}(A) = \operatorname{rank}([A|\mathbf{b}]) < n$. The number of free variables is $n - \operatorname{rank}(A)$.

---

## 01.5 Homogeneous Linear Systems

!!! definition "Definition 01.5 (Homogeneous System)"
    A system of the form $Ax = 0$ is called a homogeneous system. It always has the **trivial solution** $x = 0$.
    - A homogeneous system has non-trivial solutions if and only if $\operatorname{rank}(A) < n$.
    - If $m < n$ (fewer equations than variables), $Ax=0$ always has a non-trivial solution.

---

## Exercises

1. **[Consistency] Determine if the system $x+y=1, x+y=2$ has a solution.**
   ??? success "Solution"
       No solution (Inconsistent). The augmented matrix is $\begin{pmatrix} 1 & 1 & | & 1 \\ 1 & 1 & | & 2 \end{pmatrix}$. Row reduction gives $\begin{pmatrix} 1 & 1 & | & 1 \\ 0 & 0 & | & 1 \end{pmatrix}$. Since $\operatorname{rank}(A)=1 \neq \operatorname{rank}(\tilde{A})=2$, there is no solution.

2. **[Elimination] Use Gaussian elimination to solve: $x+y=3, x-y=1$.**
   ??? success "Solution"
       Adding the equations gives $2x=4 \implies x=2$. Substitution gives $2+y=3 \implies y=1$. The solution is $(2, 1)$.

3. **[RREF] Transform $\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$ into RREF.**
   ??? success "Solution"
       $R_2 - 2R_1 \to R_2$ gives $\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix}$.

4. **[Free Variables] If a system with 3 equations and 5 variables is full rank, how many free variables are there?**
   ??? success "Solution"
       Number of free variables $= n - \operatorname{rank}(A) = 5 - 3 = 2$.

5. **[Homogeneous] Describe the solution space of $x_1 + x_2 + x_3 = 0$.**
   ??? success "Solution"
       This is a plane through the origin in $\mathbb{R}^3$. Let $x_2=s, x_3=t$, then $x_1 = -s-t$. The solution is $s(-1, 1, 0)^T + t(-1, 0, 1)^T$.

6. **[Rank] Prove: If $Ax=b$ has a solution for every $b$, $A$ must have full row rank.**
   ??? success "Solution"
       A solution existing for every $b$ means $\operatorname{Im}(A) = \mathbb{R}^m$. Thus $\operatorname{rank}(A) = m$.

7. **[Geometry] What does a contradictory row $[0, 0 | 1]$ represent geometrically?**
   ??? success "Solution"
       It represents parallel hyperplanes that do not intersect.

8. **[Parametric] Find the general solution to $x+y+z=1$.**
   ??? success "Solution"
       Let $y=s, z=t$, then $x = 1-s-t$. The general solution is $(1, 0, 0)^T + s(-1, 1, 0)^T + t(-1, 0, 1)^T$.

9. **[Structure] What is the structure of the general solution to a non-homogeneous system?**
   ??? success "Solution"
       General Solution = Particular Solution + Homogeneous Solution.

10. **[Numerical] Why is "pivoting" used in Gaussian elimination for large-scale computing?**
    ??? success "Solution"
        To minimize the magnification of rounding errors caused by dividing by very small numbers, thereby improving numerical stability.

## Chapter Summary

This chapter bridges intuitive equations and algorithmic matrices:

1.  **Algorithmic Engine**: Gaussian elimination is the universal tool for linear problems; its RREF form reveals the intrinsic structure of the system.
2.  **Structural Insight**: The existence and uniqueness of solutions are entirely determined by the comparison of the ranks of the coefficient and augmented matrices.
3.  **Spatial Foundation**: The solution set of a homogeneous system forms a subspace, while the solution set of a non-homogeneous system is a translation of that subspace (an affine space).
