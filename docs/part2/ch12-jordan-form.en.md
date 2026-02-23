# Chapter 12: The Jordan Canonical Form

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Minimal Polynomial (Ch00) · Eigenspaces · Defective Matrices

**Chapter Outline**: Generalized Eigenvectors → Jordan Blocks → Structure of the Jordan Form → Uniqueness and Existence → Computing the Jordan Form → Chain of Generalized Eigenvectors → Power of Jordan Blocks → Relation to Minimal and Characteristic Polynomials

**Extension**: The Jordan form is the closest a non-diagonalizable matrix can get to being diagonal; it reveals the complete "cycle structure" of a linear operator.

</div>

Not every matrix can be diagonalized. Matrices whose algebraic multiplicity exceeds their geometric multiplicity—**defective matrices**—cannot be represented as $PDP^{-1}$. The **Jordan Canonical Form** provides the most simplified representation possible for *any* square matrix over an algebraically closed field. By grouping **generalized eigenvectors** into chains, we form **Jordan blocks**, which capture the internal "dragging" effect of non-diagonalizable operators. This chapter explores the construction of these blocks and their role in understanding matrix stability and the behavior of linear differential equations.

---

## 12.1 Generalized Eigenvectors and Jordan Blocks

!!! definition "Definition 12.1 (Generalized Eigenvector)"
    A vector $v \neq 0$ is a generalized eigenvector of $A$ for eigenvalue $\lambda$ if $(A - \lambda I)^k v = 0$ for some $k \ge 1$.

!!! definition "Definition 12.2 (Jordan Block)"
    A Jordan block $J_k(\lambda)$ is a $k \times k$ matrix with $\lambda$ on the diagonal, 1s on the super-diagonal, and 0s elsewhere:
    $$J_k(\lambda) = \begin{pmatrix} \lambda & 1 & 0 & \dots \\ 0 & \lambda & 1 & \dots \\ \vdots & \vdots & \ddots & 1 \\ 0 & 0 & 0 & \lambda \end{pmatrix}$$

---

## Exercises

1. **[Fundamentals] Write the Jordan Form of $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$.**
   ??? success "Solution"
       The matrix is already a single $2 \times 2$ Jordan block: $J_2(2)$. It has one eigenvalue 2 and only one linearly independent eigenvector $(1, 0)^T$.

2. **[Cycle Structure] If $\lambda$ has algebraic multiplicity 3 and geometric multiplicity 1, what is the Jordan structure?**
   ??? success "Solution"
       There is exactly one Jordan block for $\lambda$, and its size must be 3: $J_3(\lambda)$.

3. **[Minimal Polynomial] How does the minimal polynomial $m(x)$ relate to the Jordan blocks?**
   ??? success "Solution"
       The degree of $(x-\lambda)$ in $m(x)$ is the size of the **largest** Jordan block associated with $\lambda$.

4. **[Calculation] Find the Jordan form of $A = \begin{pmatrix} 3 & 0 & 0 \\ 0 & 2 & 1 \\ 0 & 0 & 2 \end{pmatrix}$.**
   ??? success "Solution"
       The matrix is block diagonal. The first block is $J_1(3)$ and the second is $J_2(2)$. Thus $J = \operatorname{diag}(J_1(3), J_2(2))$.

5. **[Powers] Calculate $J_2(0)^2$.**
   ??? success "Solution"
       $J_2(0) = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$, so $J_2(0)^2 = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$. Powers of nilpotent Jordan blocks eventually vanish.

6. **[Uniqueness] Is the Jordan form of a matrix unique?**
   ??? success "Solution"
       Yes, up to the ordering of the Jordan blocks along the diagonal.

7. **[Chain] What is a "Jordan chain"?**
   ??? success "Solution"
       A sequence of vectors $\{v_1, \dots, v_k\}$ such that $Av_1 = \lambda v_1$, $Av_2 = \lambda v_2 + v_1$, etc. These vectors form the columns of the transformation matrix $P$.

8. **[Algebraic Multiplicity] Relate the algebraic multiplicity to the sum of Jordan block sizes.**
   ??? success "Solution"
       The algebraic multiplicity of $\lambda$ is the sum of the sizes of all Jordan blocks associated with $\lambda$.

9. **[Geometric Multiplicity] Relate the geometric multiplicity to the number of Jordan blocks.**
   ??? success "Solution"
       The geometric multiplicity of $\lambda$ is exactly equal to the number of Jordan blocks associated with that eigenvalue.

10. **[Functions] Why is the Jordan form important for computing $e^{At}$?**
    ??? success "Solution"
        Because $e^{At}$ of a Jordan block has a simple explicit formula involving powers of $t$. By decomposing $A$ into Jordan form, we can calculate the exponential of any matrix analytically.

## Chapter Summary

This chapter establishes the definitive canonical form for general square operators:

1. **Cycle Decomposition**: Developed generalized eigenvectors to capture the "dragging" effects of defective matrices.
2. **Structural Block**: Defined the Jordan block as the irreducible building block of non-diagonalizable systems.
3. **Polynomial Link**: Established the connection between Jordan block sizes and the multiplicities in the minimal and characteristic polynomials.
4. **Analytic Platform**: Positioned the Jordan form as the essential tool for exact analytical solutions to linear differential equations.
