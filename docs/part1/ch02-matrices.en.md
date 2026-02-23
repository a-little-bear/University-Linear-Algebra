# Chapter 02: Matrices and Matrix Operations

<div class="context-flow" markdown>

**Prerequisites**: Basic Algebraic Operations

**Chapter Outline**: Matrix Definition and Representation → Addition and Scalar Multiplication → Matrix Multiplication (Non-commutativity) → Transpose → Special Matrices (Identity, Diagonal, Symmetric) → Matrix Inverses $A^{-1}$ → Matrix Powers and Polynomials

**Extension**: Matrices are not just tables of data, but mapping operators of linear spaces.

</div>

Matrices are the computational core of linear algebra. They condense linear transformations of vectors into rectangular arrays. The introduction of matrix multiplication is more than a stack of rules; it reflects the composition of linear maps.

---

## 02.1 Matrix Operation Rules

!!! definition "Definition 02.1 (Matrix Multiplication)"
    If $A$ is an $m \times n$ matrix and $B$ is an $n \times p$ matrix, their product $C = AB$ is an $m \times p$ matrix with entries:
    $$c_{ij} = \sum_{k=1}^n a_{ik} b_{kj}$$
    Note: In general, $AB \neq BA$.

---

## Exercises

1. **[Basic Operations] Given $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$ and $B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$. Compute $AB$ and $BA$.**
   ??? success "Solution"
       $AB = \begin{pmatrix} 1(0)+2(1) & 1(1)+2(0) \\ 3(0)+4(1) & 3(1)+4(0) \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}$.
       $BA = \begin{pmatrix} 0(1)+1(3) & 0(2)+1(4) \\ 1(1)+0(3) & 1(2)+0(4) \end{pmatrix} = \begin{pmatrix} 3 & 4 \\ 1 & 2 \end{pmatrix}$.
       Thus $AB \neq BA$.

2. **[Identity] Prove that for any $n \times n$ matrix $A$, $AI = IA = A$.**
   ??? success "Solution"
       The entries $\delta_{ij}$ of the identity matrix $I$ are 1 if $i=j$ and 0 otherwise. Substituting into the multiplication formula: $(AI)_{ij} = \sum a_{ik} \delta_{kj} = a_{ij}$.

3. **[Transpose] Known that $(AB)^T = B^T A^T$. Use this to find $(A^T B)^T$.**
   ??? success "Solution"
       $(A^T B)^T = B^T (A^T)^T = B^T A$.

4. **[Symmetric] If $A$ is a symmetric matrix, prove that $A^2$ is also symmetric.**
   ??? success "Solution"
       A symmetric matrix satisfies $A^T = A$. Then $(A^2)^T = (AA)^T = A^T A^T = AA = A^2$. Thus $A^2$ is symmetric.

5. **[Inverses Intro] Find the inverse of $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$.**
   ??? success "Solution"
       Using the $2 \times 2$ formula: $A^{-1} = \frac{1}{ad-bc} \begin{pmatrix} d & -b \\ -c & a \end{pmatrix}$.
       $\det A = 1(4)-2(3) = -2$.
       $A^{-1} = \frac{1}{-2} \begin{pmatrix} 4 & -2 \\ -3 & 1 \end{pmatrix} = \begin{pmatrix} -2 & 1 \\ 1.5 & -0.5 \end{pmatrix}$.

6. **[Powers] Calculate $A^k$ where $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       $A^2 = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}, A^3 = \begin{pmatrix} 1 & 3 \\ 0 & 1 \end{pmatrix}$.
       By induction, $A^k = \begin{pmatrix} 1 & k \\ 0 & 1 \end{pmatrix}$.

7. **[Matrix Equations] If $AX = B$ and $A$ is invertible, solve for $X$.**
   ??? success "Solution"
       Left-multiply by $A^{-1}$: $A^{-1}AX = A^{-1}B \implies X = A^{-1}B$. Note that you cannot right-multiply.

8. **[Rank Intro] What is the rank of the matrix $\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$?**
   ??? success "Solution"
       The rank is 1. It has only one non-zero row (or one linearly independent column).

9. **[Orthogonality Intro] Verify that $\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$ is orthogonal (i.e., $AA^T = I$).**
   ??? success "Solution"
       $A^T = \begin{pmatrix} \cos & \sin \\ -\sin & \cos \end{pmatrix}$. Multiplication yields 1 on the diagonal ($\cos^2+\sin^2=1$) and 0 off-diagonal ($\cos\sin-\sin\cos=0$). Thus $AA^T = I$.

10. **[Diagonal] Is the product of two diagonal matrices still diagonal? What is the rule for its entries?**
    ??? success "Solution"
        Yes. The result is a diagonal matrix where entries are multiplied term-wise: $\operatorname{diag}(a_i) \operatorname{diag}(b_i) = \operatorname{diag}(a_i b_i)$.

## Chapter Summary

Matrix operations build the computational syntax of linear algebra:

1. **Non-commutativity**: This is the most fundamental difference between matrix and scalar multiplication.
2. **Structural Preservation**: Transpose and inverse operations maintain the internal algebraic logic.
3. **Operational Simplification**: Special matrices (identity, diagonal) greatly simplify the analysis of complex systems.
