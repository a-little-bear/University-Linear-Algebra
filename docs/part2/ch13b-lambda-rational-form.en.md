# Chapter 13B: λ-Matrices and Rational Canonical Form

<div class="context-flow" markdown>

**Prerequisites**: Polynomial Algebra (Ch0) · Jordan Form (Ch12) · Smith Normal Form (Ch48)

**Chapter Outline**: Definition of λ-Matrices → Determinantal Divisors → Invariant Factors → Elementary Divisors → Necessary and Sufficient Conditions for Matrix Similarity → Frobenius Rational Canonical Form → Definition and Properties of Minimal Polynomials → Relation to Diagonalization

**Extension**: Rational Canonical Form solves the classification problem of matrix similarity over general fields (not just complex numbers), generalizing the Jordan form.

</div>

In the study of Jordan form, we rely on the algebraic closure of the complex field. However, over the real field or general fields, characteristic polynomials may not split completely. **λ-matrix** theory, by studying the Smith Normal Form of polynomial matrices, reveals invariants of similarity transformations under any field.

---

## 13B.1 Core Concepts

!!! definition "Definition 13B.1 (Invariant Factors)"
    For the characteristic matrix $\lambda I - A$, the non-zero polynomials $d_1(\lambda), \dots, d_n(\lambda)$ on the diagonal of its Smith Normal Form are called the **invariant factors** of $A$. They satisfy $d_i \mid d_{i+1}$.

!!! theorem "Theorem 13B.3 (Similarity Criterion)"
    Two matrices $A, B$ are similar if and only if their characteristic matrices $\lambda I - A$ and $\lambda I - B$ have the same invariant factors (or elementary divisors).

---

## Exercises

1. **[Fundamentals] Write the characteristic matrix $\lambda I - A$ for $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$.**
   ??? success "Solution"
       $\lambda I - A = \begin{pmatrix} \lambda & -1 \\ 0 & \lambda \end{pmatrix}$.

2. **[Invariant Factors Calculation] Find the invariant factors of $A = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$.**
   ??? success "Solution"
       Characteristic matrix: $\begin{pmatrix} \lambda-2 & 0 \\ 0 & \lambda-2 \end{pmatrix}$.
       GCD of 1st-order minors: $D_1 = \lambda-2$.
       GCD of 2nd-order minors (determinant): $D_2 = (\lambda-2)^2$.
       Invariant factors: $d_1 = D_1 = \lambda-2$ and $d_2 = D_2/D_1 = \lambda-2$.

3. **[Minimal Polynomial] How does the minimal polynomial $m(\lambda)$ relate to invariant factors?**
   ??? success "Solution"
       The minimal polynomial equals the **last invariant factor** $d_n(\lambda)$ in the Smith form. It contains all eigenvalues as roots, and its power corresponds to the largest Jordan block size.

4. **[Diagonalization Criterion] State the condition for diagonalizability using the minimal polynomial.**
   ??? success "Solution"
       A matrix is diagonalizable if and only if its minimal polynomial $m(\lambda)$ has no multiple roots (splits into distinct linear factors).

5. **[Elementary Divisors] Given invariant factors $d_1=1, d_2=(\lambda-1), d_3=(\lambda-1)^2(\lambda-2)$. Find the elementary divisors.**
   ??? success "Solution"
       Elementary divisors are the highest power terms in the irreducible factorization of the invariant factors.
       In this case: $(\lambda-1), (\lambda-1)^2, (\lambda-2)$.

6. **[Calculation] Write the companion matrix $C(p)$ for $p(\lambda) = \lambda^2 + a\lambda + b$.**
   ??? success "Solution"
       $C(p) = \begin{pmatrix} 0 & -b \\ 1 & -a \end{pmatrix}$. Its characteristic polynomial is exactly $p(\lambda)$.

7. **[Rational Form] What is the main difference between Rational Canonical Form and Jordan Form?**
   ??? success "Solution"
       Rational form does not require finding roots over the field and is composed of companion blocks; Jordan form requires the polynomial to split. Rational form exists and is unique over *any* field.

8. **[Cayley-Hamilton] Use invariant factors to prove the Cayley-Hamilton Theorem.**
   ??? success "Solution"
       The characteristic poly is $f(\lambda) = d_1 d_2 \dots d_n$, and the minimal poly is $m(\lambda) = d_n$. Clearly $m(\lambda) \mid f(\lambda)$. Since $m(A)=O$, it follows that $f(A)=O$.

9. **[Matrix Rank] If the Smith form of $\lambda I - A$ has $r$ non-zero diagonal entries, what does it mean?**
   ??? success "Solution"
       For the characteristic matrix of an $n \times n$ matrix, the determinant is a non-zero polynomial, so the rank is always $n$ (over the field of fractions). All $n$ diagonal entries are non-zero polynomials.

10. **[Application] Why are invariant factors better than eigenvalues for studying the similarity of real matrices?**
    ??? success "Solution"
        Because invariant factors are polynomials that can be computed via the Euclidean algorithm within $\mathbb{R}$, without dealing with imaginary eigenvalues. Two real matrices are similar over $\mathbb{C}$ iff they are similar over $\mathbb{R}$.

## Chapter Summary

λ-matrix theory is the pinnacle of the structural theory of linear algebra:

1. **Unified View**: Unifies all similarity invariants (rank, eigenvalues, Jordan blocks) via Smith Normal Form.
2. **Field Independence**: Rational Canonical Form provides a classification scheme independent of field extensions.
3. **Minimal Polynomial**: Establishes the final link between matrix nilpotency and the algebraic properties of polynomials.
