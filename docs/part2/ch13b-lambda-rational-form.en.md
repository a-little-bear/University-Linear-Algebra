# Chapter 13B: λ-matrices and Rational Canonical Form

<div class="context-flow" markdown>

**Prerequisites**: Polynomial Algebra (Ch00) · Matrix Algebra (Ch02) · Jordan Canonical Form (Ch12)

**Chapter Outline**: Definition of $\lambda$-matrices → Elementary Operations and Equivalence → Smith Normal Form → Invariant Factors → Elementary Divisors → Necessary and Sufficient Conditions for Similarity → Companion Matrices → Rational Canonical Form (RCF) → Structural Comparison (Jordan vs. Rational)

**Extension**: The Rational Canonical Form solves the problem of matrix canonical forms over general fields (e.g., $\mathbb{Q}$) without requiring the field to be algebraically closed (e.g., $\mathbb{C}$); it is a concrete application of the structure theorem for modules over a PID.

</div>

While the Jordan Canonical Form is theoretically perfect, its construction depends on the existence of eigenvalues within the underlying field. Working over the field of rational numbers $\mathbb{Q}$, the characteristic polynomial may not factor into linear terms. The Rational Canonical Form (RCF) overcomes this limitation by providing a universal standard form applicable over any field through the decomposition of polynomial rings.

---

## 13B.1 λ-matrices and Smith Normal Form

!!! definition "Definition 13B.1 (λ-matrix)"
    A matrix whose entries are polynomials in $\lambda$ is called a **$\lambda$-matrix** (or a matrix polynomial).

!!! theorem "Theorem 13B.1 (Smith Normal Form)"
    Every $n \times n$ $\lambda$-matrix $A(\lambda)$ can be transformed via elementary operations into a unique diagonal form:
    $$S(\lambda) = \operatorname{diag}(d_1(\lambda), d_2(\lambda), \ldots, d_r(\lambda), 0, \ldots, 0)$$
    where each $d_i(\lambda)$ is a monic polynomial such that $d_i(\lambda) \mid d_{i+1}(\lambda)$. These polynomials are called the **invariant factors** of $A(\lambda)$.

---

## 13B.2 Conditions for Similarity

!!! theorem "Theorem 13B.2 (Similarity Theorem)"
    Two $n \times n$ matrices $A$ and $B$ are similar if and only if their characteristic matrices $\lambda I - A$ and $\lambda I - B$ have the same Smith Normal Form (i.e., identical invariant factors).

---

## 13B.3 Companion Matrices and Rational Canonical Form

!!! definition "Definition 13B.2 (Companion Matrix $C(p)$)"
    For a monic polynomial $p(\lambda) = \lambda^k + a_{k-1}\lambda^{k-1} + \cdots + a_0$, the **companion matrix** is:
    $$C(p) = \begin{pmatrix} 0 & 0 & \cdots & -a_0 \\ 1 & 0 & \cdots & -a_1 \\ 0 & 1 & \cdots & -a_2 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & -a_{k-1} \end{pmatrix}$$
    **Property**: The characteristic and minimal polynomials of $C(p)$ are both $p(\lambda)$.

!!! theorem "Theorem 13B.3 (Rational Canonical Form)"
    Every square matrix $A$ is similar to a block diagonal matrix whose blocks are the companion matrices of its invariant factors:
    $$R = \operatorname{diag}(C(d_1), C(d_2), \ldots, C(d_k))$$
    This is known as the **Rational Canonical Form** of $A$.

---

## Exercises

****

??? success "Solution"
    
       Thus $d_1 = D_1 = 1$ and $d_2 = D_2/D_1 = \lambda^2$. The invariant factors are $1, \lambda^2$.

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

## Chapter Summary

The Rational Canonical Form provides a universal characterization of matrix similarity classes over any field:

1.  **Field Independence**: By utilizing companion blocks of irreducible polynomials, RCF eliminates the dependency on the complex field, becoming the core of operator theory in abstract algebra.
2.  **Polynomial Logic**: The theory of invariant factors reveals the deep module structure behind the characteristic matrix $\lambda I - A$, establishing the ultimate algebraic criterion for similarity.
3.  **Computational Precision**: Compared to the instability of eigenvalue calculation, the construction of Smith forms based on elementary operations provides a robust algorithmic foundation for exact algebra software.
