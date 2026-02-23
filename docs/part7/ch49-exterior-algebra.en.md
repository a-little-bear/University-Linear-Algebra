# Chapter 49: Exterior Algebra and Grassmannians

<div class="context-flow" markdown>

**Prerequisites**: Multilinear Algebra & Tensors (Ch21) · Determinants (Ch03) · Vector Spaces (Ch04)

**Chapter Outline**: Alternating Multilinear Forms → Definition and Properties of the Wedge Product ($\wedge$) → Exterior Power Space $\Lambda^k(V)$ (Dimension and Basis) → Compound Matrices ($C_k(A)$) → The Essence of the Binet-Cauchy Formula → Definition of the Grassmannian $Gr(k, V)$ → Plücker Embedding and Plücker Coordinates → Exterior Algebra as the Theoretical Home of Determinants → Applications: Differential Forms and Combinatorial Geometry

**Extension**: Exterior algebra provides the foundation for modern Differential Geometry (Differential Forms) and Topology; Plücker coordinates transform the study of subspaces into the study of algebraic varieties in projective space, serving as a classic example in Algebraic Geometry.

</div>

What is the algebraic essence of the determinant? Why is it alternating? **Exterior Algebra** (also known as Grassmann Algebra) provides the ultimate answer. By introducing the **Wedge Product** ($\wedge$), we can formally describe "oriented volume" and treat a $k$-dimensional subspace as a single algebraic object. This chapter reveals the geometric-algebraic structures hidden behind the determinant.

---

## 49.1 The Wedge Product and Exterior Powers

!!! definition "Definition 49.1 (The Wedge Product $\wedge$)"
    Let $V$ be a vector space. For $u, v \in V$, the **wedge product** $u \wedge v$ satisfies:
    1.  **Bilinearity**: Distributive over addition and homogeneous under scalar multiplication.
    2.  **Antisymmetry**: $v \wedge u = -(u \wedge v)$.
    3.  **Nilpotency**: $v \wedge v = 0$.

!!! theorem "Theorem 49.1 (Dimension of Exterior Powers)"
    The space spanned by the wedge products of $k$ vectors from $V$ is the **$k$-th exterior power**, denoted $\Lambda^k(V)$.
    If $\dim V = n$, then $\dim \Lambda^k(V) = \binom{n}{k}$.
    When $k=n$, $\Lambda^n(V)$ is 1-dimensional, which is the essence of the determinant.

---

## 49.2 Compound Matrices and Binet-Cauchy

!!! definition "Definition 49.2 (Compound Matrix $C_k(A)$)"
    For an $m \times n$ matrix $A$, the **$k$-th compound matrix** $C_k(A)$ is a $\binom{m}{k} \times \binom{n}{k}$ matrix whose entries are all possible $k \times k$ minors of $A$.

!!! theorem "Theorem 49.2 (Operator Perspective)"
    The compound matrix $C_k(A)$ is the representation of the linear transformation induced by $A$ on the exterior space $\Lambda^k(V)$.
    This leads to the **Binet-Cauchy Formula**: $C_k(AB) = C_k(A)C_k(B)$.

---

## 49.3 The Grassmannian and Plücker Embedding

!!! definition "Definition 49.3 (Grassmannian)"
    The **Grassmannian** $Gr(k, V)$ is the set of all $k$-dimensional subspaces of $V$. It is a compact smooth manifold.

!!! technique "The Plücker Embedding"
    A subspace $W = \operatorname{span}\{w_1, \ldots, w_k\}$ is mapped to a point in the exterior power space:
    $$\Phi(W) = [w_1 \wedge w_2 \wedge \cdots \wedge w_k] \in \mathbb{P}(\Lambda^k(V))$$
    The resulting coordinates are called **Plücker coordinates**. They must satisfy a set of quadratic relations known as the Plücker relations.

---

## Exercises

1.  **[Basics] In $\mathbb{R}^3$, calculate $(e_1 + e_2) \wedge (e_2 + e_3)$.**
    ??? success "Solution"
        $= e_1 \wedge e_2 + e_1 \wedge e_3 + e_2 \wedge e_2 + e_2 \wedge e_3 = e_1 \wedge e_2 - e_3 \wedge e_1 + e_2 \wedge e_3$.

2.  **[Dimension] If $\dim V = 4$, what is the dimension of $\Lambda^2(V)$?**
    ??? success "Solution"
        $\binom{4}{2} = 6$.

3.  **[Independence] Prove that $v_1, \ldots, v_k$ are linearly independent iff $v_1 \wedge \cdots \wedge v_k \neq 0$.**
    ??? success "Solution"
        If dependent, one vector is a combination of others; by antisymmetry, the product is 0. If independent, they form a part of a basis, and their product is a basis vector of $\Lambda^k(V)$, hence non-zero.

4.  **[Determinant] Prove that for an $n \times n$ matrix $A$, $Av_1 \wedge \cdots \wedge Av_n = \det(A)(v_1 \wedge \cdots \wedge v_n)$.**
    ??? success "Solution"
        This is a property of the 1-dimensional space $\Lambda^n(V)$. Any linear operator acts as a scalar multiplier on this space, and that scalar is the determinant.

5.  **[Compound] Write the 2nd compound matrix of $\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$.**
    ??? success "Solution"
        It has only one element, which is the $2 \times 2$ minor (determinant): $(-2)$.

6.  **[Plücker] Find the Plücker coordinates of $(1, 0, 0, 0) \wedge (0, 1, 0, 0)$ in $\mathbb{P}(\Lambda^2(\mathbb{R}^4))$.**
    ??? success "Solution"
        In the basis $\{e_i \wedge e_j\}$, only the $e_1 \wedge e_2$ component is 1. The coordinates are $(1, 0, 0, 0, 0, 0)$.

7.  **[Associativity] Prove $(u \wedge v) \wedge w = u \wedge (v \wedge w)$.**
    ??? success "Solution"
        The wedge product is associative, guaranteed by the properties of the tensor product and the antisymmetrization operator.

8.  **[Trace] Prove $\operatorname{tr}(C_k(A))$ is the sum of all possible products of $k$ eigenvalues of $A$.**
    ??? success "Solution"
        The eigenvalues of $C_k(A)$ are exactly the products of $k$ eigenvalues of $A$. The trace is the sum of these products (the elementary symmetric polynomials in the eigenvalues).

9.  **[Geometry] Why is exterior algebra the foundation for differential forms?**
    ??? success "Solution"
        Differential forms are sections of the exterior power of the cotangent bundle. The wedge product corresponds to the multiplication of these forms.

****

??? success "Solution"
    

## Chapter Summary

Exterior algebra is the ultimate language for determinants and geometric subspaces:

1.  **Power of Antisymmetry**: By introducing the $\wedge$ operator, exterior algebra precisely translates geometric "orientation" and "volume" into algebraic sign-switching, establishing a standard framework for oriented geometric quantities.
2.  **Dimension Elevation**: Compound matrix theory proves that the action of linear transformations on higher-order tensor spaces maintains perfect structure (Binet-Cauchy), serving as a powerful tool for studying complex determinant identities.
3.  **Subspaces as Objects**: The Grassmannian and Plücker embedding transform dynamic subspace selection into static projective points, marking the transition from elementary linear algebra to algebraic geometry and manifold theory.
