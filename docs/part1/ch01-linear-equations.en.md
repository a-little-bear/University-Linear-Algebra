# Chapter 01: Systems of Linear Equations

<div class="context-flow" markdown>

**Prerequisites**: Basic Algebra

**Chapter Outline**: Definition of Linear Systems → Augmented Matrix → Elementary Row Operations → Echelon Forms (REF/RREF) → Gaussian Elimination → Consistency and Uniqueness of Solutions

**Extension**: Gaussian elimination is not just a tool for solving equations, but also the basis for calculating matrix rank and inverses.

</div>

Systems of linear equations are the starting point of linear algebra. Almost all engineering, physical, and economic models eventually boil down to solving $Ax = b$. Through elementary operations, we can simplify complex systems into echelon forms that reveal their solution structure.

---

## 01.1 Gaussian Elimination

!!! definition "Definition 01.1 (Elementary Row Operations)"
    1. Swap two rows;
    2. Multiply a row by a non-zero constant;
    3. Add a multiple of one row to another.
    These operations do not change the solution set of the system.

---

## Exercises

1. **[Fundamentals] Determine if the system $x+y=1, x+y=2$ has a solution. Explain.**
   ??? success "Solution"
       No solution (inconsistent). The left sides are identical but the right sides are different, representing two parallel lines that never intersect.

2. **[Elimination] Use Gaussian elimination to solve: $x+y=3, x-y=1$.**
   ??? success "Solution"
       Add them: $2x=4 \implies x=2$. Substitute: $2+y=3 \implies y=1$. The solution is $(2, 1)$.

3. **[Matrix Form] Write the augmented matrix for the system $2x-y=5, x+3y=0$.**
   ??? success "Solution"
       $\begin{pmatrix} 2 & -1 & | & 5 \\ 1 & 3 & | & 0 \end{pmatrix}$.

4. **[Consistency] If the echelon form of an augmented matrix contains a row like $[0, 0, 0 | 1]$, what does it mean?**
   ??? success "Solution"
       It means the system is inconsistent (no solution). That row represents $0x + 0y + 0z = 1$, which is logically impossible.

5. **[Free Variables] In RREF form, what do columns without pivots correspond to?**
   ??? success "Solution"
       They correspond to **free variables**. These variables can take any value (usually denoted by parameters $t, s$), leading to infinitely many solutions.

6. **[Homogeneous Systems] Prove: A homogeneous system $Ax = 0$ always has at least one solution.**
   ??? success "Solution"
       Yes, it always has the trivial solution $x=0$. If the number of variables exceeds the number of equations, it must have non-trivial solutions.

7. **[Calculation] Transform $\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$ into Reduced Row Echelon Form (RREF).**
   ??? success "Solution"
       Subtract 2 times row 1 from row 2 to get $\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix}$. This is the RREF.

8. **[Application] Are the equations derived from Kirchhoff's laws for a circuit typically linear?**
   ??? success "Solution"
       Yes. Kirchhoff's voltage and current laws involve linear superpositions of voltage drops and currents, so the resulting systems are linear.

9. **[Parametric Solution] Solve $x+y+z=1$.**
   ??? success "Solution"
       There are two free variables. Let $y=s, z=t$, then $x = 1-s-t$. The general solution is $\begin{pmatrix} x \\ y \\ z \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix} + s\begin{pmatrix} -1 \\ 1 \\ 0 \end{pmatrix} + t\begin{pmatrix} -1 \\ 0 \\ 1 \end{pmatrix}$.

10. **[Uniqueness] What is the condition for a unique solution (in terms of pivots)?**
    ??? success "Solution"
        In the echelon form, every column (except the last constant column) contains a pivot, and there are no inconsistent rows.

## Chapter Summary

This chapter makes the jump from intuitive equations to abstract matrices:

1. **Algorithm Core**: Gaussian elimination is the universal algorithm for solving linear problems.
2. **Structural Insight**: Through echelon forms, we can clearly distinguish between no solution, a unique solution, and infinitely many solutions.
3. **Parametrization**: Introduced free variables to describe the geometry of high-dimensional solution spaces.
