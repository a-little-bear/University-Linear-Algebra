# Chapter 02: Matrices and Matrix Algebra

<div class="context-flow" markdown>

**Prerequisites**: Linear Equations (Ch1) · Vector Spaces (Ch4) · Function Composition

**Chapter Outline**: Definition of Matrices → Matrix Operations (Addition, Scaling, Multiplication) → Properties of Matrix Multiplication (Non-commutativity) → Special Matrices (Identity, Zero, Symmetric, Diagonal) → Matrix Inversion → Transpose and Adjoint → Partitioned Matrices → Matrix Algebra Laws

**Extension**: Matrices are the concrete representations of linear operators; matrix multiplication is the algebraic expression of the composition of linear transformations.

</div>

Matrices provide a compact and powerful notation for handling linear data. While they appear as simple rectangular arrays of numbers, they follow a rigorous algebraic structure. Matrix multiplication, in particular, is defined to match the composition of linear maps, leading to its non-commutative nature ($AB \neq BA$). This chapter transitions from solving equations to manipulating operators, establishing the rules for matrix inversion and the properties of the transpose.

---

## 02.1 Matrix Multiplication and Inverses

!!! definition "Definition 02.1 (Matrix Multiplication)"
    Let $A \in M_{m \times n}$ and $B \in M_{n \times p}$. The product $C = AB$ is an $m \times p$ matrix where:
    $$c_{ij} = \sum_{k=1}^n a_{ik} b_{kj}$$
    This row-by-column operation corresponds to applying operator $B$ followed by operator $A$.

!!! theorem "Theorem 02.1 (Properties of the Inverse)"
    If $A$ and $B$ are invertible, then $(AB)^{-1} = B^{-1} A^{-1}$. The order of factors is reversed, mirroring the "socks-then-shoes" logic of undoing compositions.

---

## Exercises

1. **[Fundamentals] Compute $AB$ for $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$ and $B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.**
   ??? success "Solution"
       $AB = \begin{pmatrix} 1(0)+2(1) & 1(1)+2(0) \\ 3(0)+4(1) & 3(1)+4(0) \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}$. Multiplication by $B$ effectively swapped the columns of $A$.

2. **[Non-commutativity] Provide a $2 \times 2$ example where $AB \neq BA$.**
   ??? success "Solution"
       Let $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}, B = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}$. $AB = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$, but $BA = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$. Matrix algebra is generally non-commutative.

3. **[Inversion] Calculate the inverse of $\begin{pmatrix} a & b \\ c & d \end{pmatrix}$ using the formula.**
   ??? success "Solution"
       $A^{-1} = \frac{1}{ad-bc} \begin{pmatrix} d & -b \\ -c & a \end{pmatrix}$, provided $ad-bc \neq 0$.

4. **[Transpose] Prove that $(AB)^T = B^T A^T$.**
   ??? success "Solution"
       The $(i,j)$ entry of $(AB)^T$ is the $(j,i)$ entry of $AB$, which is $\sum_k a_{jk} b_{ki}$. The $(i,j)$ entry of $B^T A^T$ is $\sum_k (B^T)_{ik} (A^T)_{kj} = \sum_k b_{ki} a_{jk}$. The expressions are identical.

5. **[Block Matrices] Multiply $\begin{pmatrix} I & A \\ 0 & I \end{pmatrix} \begin{pmatrix} I & -A \\ 0 & I \end{pmatrix}$ using block matrix rules.**
   ??? success "Solution"
       The product is $\begin{pmatrix} I \cdot I + A \cdot 0 & I(-A) + A \cdot I \\ 0 \cdot I + I \cdot 0 & 0(-A) + I \cdot I \end{pmatrix} = \begin{pmatrix} I & 0 \\ 0 & I \end{pmatrix} = I$. This shows the inverse of a shear-type block matrix.

6. **[Trace] Show that $\operatorname{tr}(AB) = \operatorname{tr}(BA)$ for any matrices where the products are defined.**
   ??? success "Solution"
       $\operatorname{tr}(AB) = \sum_i (AB)_{ii} = \sum_i \sum_j a_{ij} b_{ji}$. $\operatorname{tr}(BA) = \sum_j (BA)_{jj} = \sum_j \sum_i b_{ji} a_{ij}$. Both are sums over the same $n \times m$ set of products $a_{ij} b_{ji}$.

7. **[Orthogonality] Define an orthogonal matrix and state its property regarding the inverse.**
   ??? success "Solution"
       A matrix $Q$ is orthogonal if $Q^T Q = I$. This implies $Q^{-1} = Q^T$. These matrices represent rotations and reflections.

8. **[Nilpotency] If $A^2 = 0$, what is $(I-A)^{-1}$?**
   ??? success "Solution"
       $(I-A)(I+A) = I - A^2 = I$. Thus $(I-A)^{-1} = I+A$. This is a truncated geometric series.

9. **[Symmetry] If $A$ and $B$ are symmetric, is $AB$ symmetric?**
   ??? success "Solution"
       $(AB)^T = B^T A^T = BA$. Thus $AB$ is symmetric if and only if $A$ and $B$ commute ($AB = BA$).

10. **[Rank] What is the rank of the product $AB$ compared to the ranks of $A$ and $B$?**
    ??? success "Solution"
        $\operatorname{rank}(AB) \le \min(\operatorname{rank}(A), \operatorname{rank}(B))$. Multiplication can only collapse dimensions, never create them.

## Chapter Summary

This chapter defines the calculus of linear operators:

1. **Operational Logic**: Formalized matrix multiplication as the algebraic counterpart to function composition.
2. **Inversion Criteria**: Established the necessary and sufficient conditions for undoing linear transformations.
3. **Symmetry and Transposition**: Analyzed the role of the adjoint and transpose in preserving and reflecting geometric properties.
4. **Partitioned Calculus**: Developed block matrix arithmetic for managing high-dimensional systems efficiently.
