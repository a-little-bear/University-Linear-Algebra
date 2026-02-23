# Chapter 48: Modules over Principal Ideal Domains

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Rings and Fields · Determinants (Ch3)

**Chapter Outline**: Rings and Modules → Principal Ideal Domains (PIDs) → Smith Normal Form (SNF) → Unimodular Equivalence → Invariant Factors → Structure Theorem for Finitely Generated Modules over PIDs → Application to Jordan Canonical Form

**Extension**: Smith Normal Form is the unifying algebraic tool for understanding both the Jordan form of matrices and the structure of finitely generated Abelian groups.

</div>

The study of matrices over fields generalizes naturally to matrices over rings. When the ring is a **Principal Ideal Domain** (PID), such as the ring of integers $\mathbb{Z}$ or the ring of polynomials $F[x]$, we can define a canonical form known as the **Smith Normal Form**. This form provides a complete set of invariants for the equivalence of matrices under unimodular transformations and is central to the classification of finitely generated modules.

---

## 48.1 Smith Normal Form (SNF)

!!! definition "Definition 48.1 (Unimodular Matrix)"
    A square matrix $U$ over a ring $R$ is **unimodular** if its determinant is a unit in $R$. This is equivalent to saying $U^{-1}$ exists and has entries in $R$.

!!! theorem "Theorem 48.1 (Existence of SNF)"
    Let $A$ be an $m 	imes n$ matrix over a PID $R$. Then there exist unimodular matrices $P \in M_m(R)$ and $Q \in M_n(R)$ such that $PAQ = D$, where $D$ is a diagonal-like matrix:
    $$D = \begin{pmatrix} d_1 & & & 0 \ & d_2 & & \ & & \ddots & \ 0 & & & d_r \end{pmatrix}, \quad d_i \mid d_{i+1}$$
    The elements $d_1, \dots, d_r$ are unique up to multiplication by units and are called the **invariant factors** of $A$.

---

## Exercises

1. **[Invariant Factors] Compute the Smith Normal Form of the integer matrix $A = \begin{pmatrix} 2 & 4 \ 4 & 8 \end{pmatrix}$ over $\mathbb{Z}$.**
   ??? success "Solution"
       The greatest common divisor of all entries is $d_1 = \gcd(2, 4, 4, 8) = 2$. The determinant is 0, so the rank is 1. The SNF is $\begin{pmatrix} 2 & 0 \ 0 & 0 \end{pmatrix}$.

2. **[Determinantal Divisors] Define the $k$-th determinantal divisor $\delta_k(A)$ and its relation to the invariant factors $d_i$.**
   ??? success "Solution"
       $\delta_k(A)$ is the GCD of all $k 	imes k$ minors of $A$. The invariant factors are given by $d_k = \delta_k(A) / \delta_{k-1}(A)$ (with $\delta_0 = 1$).

3. **[Polynomial Ring] Find the SNF of $A(x) = \begin{pmatrix} x-1 & 0 \ 0 & (x-1)^2 \end{pmatrix}$ over $\mathbb{C}[x]$.**
   ??? success "Solution"
       The matrix is already in SNF since $(x-1)$ divides $(x-1)^2$. The invariant factors are $d_1(x) = x-1$ and $d_2(x) = (x-1)^2$.

4. **[Abelian Groups] Use SNF to find the structure of the Abelian group $G = \mathbb{Z}^2 / \operatorname{Im}(A)$ where $A = \begin{pmatrix} 2 & 0 \ 0 & 3 \end{pmatrix}$.**
   ??? success "Solution"
       The SNF of $A$ is $\begin{pmatrix} 1 & 0 \ 0 & 6 \end{pmatrix}$ (since $\gcd(2,3)=1$ and $\det=6$). Thus $G \cong \mathbb{Z}/1\mathbb{Z} \oplus \mathbb{Z}/6\mathbb{Z} \cong \mathbb{Z}_6$.

5. **[Rational Form] Explain how the SNF of the characteristic matrix $xI - A$ determines the Rational Canonical Form of $A$.**
   ??? success "Solution"
       The invariant factors $d_i(x)$ of $xI - A$ are precisely the invariant factors of the matrix $A$. The Rational Canonical Form is the block diagonal matrix of the companion matrices of these $d_i(x)$.

6. **[Unimodular Equivalence] Prove that two matrices $A, B$ are unimodularly equivalent iff they have the same invariant factors.**
   ??? success "Solution"
       Necessity follows from the fact that unimodular operations preserve the GCD of minors (determinantal divisors). Sufficiency follows from the existence and uniqueness of the SNF.

7. **[Algorithm] Outline the row and column operations required to zero out entries in the first row and column of a matrix over $\mathbb{Z}$.**
   ??? success "Solution"
       Use the Euclidean algorithm. By repeated subtractions (or division steps), the entry with the smallest absolute value is moved to $(1,1)$ and used to zero out others. If a remainder exists, repeat until the GCD is at $(1,1)$.

8. **[Jordan Form] How does the SNF of $xI - A$ relate to the Jordan blocks of $A$?**
   ??? success "Solution"
       The invariant factors $d_i(x)$ can be decomposed into elementary divisors $(x-\lambda_j)^{k_{ij}}$. Each elementary divisor corresponds to a Jordan block of size $k_{ij}$ for eigenvalue $\lambda_j$.

9. **[Module Structure] State the Structure Theorem for finitely generated modules over a PID in terms of SNF.**
   ??? success "Solution"
       Any such module $M$ is isomorphic to $R^k \oplus R/(d_1) \oplus \dots \oplus R/(d_r)$, where $R^k$ is the free part and the $d_i$ are the invariant factors of the relation matrix.

10. **[Rank] Prove that the number of non-zero invariant factors in the SNF of $A$ is equal to the rank of $A$.**
    ??? success "Solution"
        Unimodular equivalence preserves the rank because $P$ and $Q$ are invertible over the quotient field. In diagonal form, the rank is clearly the number of non-zero diagonal entries.

## Chapter Summary

This chapter establishes the Smith Normal Form as the fundamental structural decomposition for matrices over PIDs:

1. **Algebraic Invariants**: Defined invariant factors and determinantal divisors as the complete set of invariants for unimodular equivalence.
2. **Canonical Decomposition**: Developed the SNF algorithm to diagonalize matrices via row and column operations in general rings.
3. **Module Theory**: Linked matrix theory to the classification of finitely generated modules, providing a unified framework for Abelian groups and linear operators.
4. **Spectral Link**: Demonstrated how the SNF of the characteristic matrix generates both the Rational and Jordan canonical forms.
