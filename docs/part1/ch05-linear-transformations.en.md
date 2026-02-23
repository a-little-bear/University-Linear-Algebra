# Chapter 05: Linear Transformations

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch4) · Matrix Operations (Ch2)

**Chapter Outline**: Definition of Linear Transformations → Kernel and Image → Matrix Representation $T \leftrightarrow A$ → Change of Basis and Similarity → Homomorphism and Isomorphism → Operator Forms of Projection, Reflection, and Rotation

**Extension**: Linear transformations provide an "active" perspective, describing how space is stretched, twisted, or compressed.

</div>

If vector spaces are the static stage, linear transformations are the movements performed on that stage. Every linear transformation, given a basis, can be uniquely represented by a matrix.

---

## 05.1 Definitions and Core Properties

!!! definition "Definition 05.1 (Linear Transformation)"
    A mapping $T: V \to W$ is linear if for all $u, v \in V$ and scalar $c$:
    1. $T(u+v) = T(u) + T(v)$
    2. $T(cv) = cT(v)$

!!! theorem "Theorem 05.3 (Change of Basis Formula)"
    Let $A$ be the matrix of $T$ in basis $B$, and $A'$ the matrix in $B'$. If $P$ is the transition matrix from $B'$ to $B$, then:
    $$A' = P^{-1} A P$$
    Such matrices are called **similar matrices**.

---

## Exercises

1. **[Criteria] Is the mapping $T(x, y) = (x+1, y)$ a linear transformation? Explain.**
   ??? success "Solution"
       No. A linear transformation must satisfy $T(0) = 0$. Here $T(0, 0) = (1, 0) \neq (0, 0)$. It also fails the scalar multiplication property.

2. **[Matrix Representation] Let $T: \mathbb{R}^2 \to \mathbb{R}^2$ be the transformation that projects vectors onto the $x$-axis. Write its matrix in the standard basis.**
   ??? success "Solution"
       $T(1, 0) = (1, 0)$ and $T(0, 1) = (0, 0)$.
       Thus the matrix is $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$.

3. **[Kernel and Image] Let $T$ have the matrix $A = \begin{pmatrix} 1 & 2 \\ 2 & 2 \end{pmatrix}$. Calculate $\dim(\ker T)$.**
   ??? success "Solution"
       $\det A = 2-4 = -2 \neq 0$. The matrix has full rank, so the kernel contains only the zero vector.
       $\dim(\ker T) = 0$.

4. **[Basis Change] Given $A = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$. If the change of basis matrix is $P = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$, find the similar matrix $A'$.**
   ??? success "Solution"
       $P^{-1} = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$.
       $A' = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 2 & -3 \\ 0 & 3 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 2 & -1 \\ 0 & 3 \end{pmatrix}$.

5. **[Differential Operator] Prove that the derivative $D: P_2 \to P_1$ ($D(p) = p'$) is a linear transformation.**
   ??? success "Solution"
       By calculus rules: $(f+g)' = f' + g'$ and $(cf)' = cf'$. Since it satisfies additivity and scalar homogeneity, $D$ is linear.

6. **[Isomorphism] Prove: Vector spaces of the same dimension are isomorphic.**
   ??? success "Solution"
       Let $V, W$ have dimension $n$. Choose bases for each and define a linear map mapping one basis to the other. This map is a bijection and preserves operations, so $V \cong W$.

7. **[Rotation Matrix] Write the matrix for a $90^\circ$ counter-clockwise rotation in $\mathbb{R}^2$.**
   ??? success "Solution"
       $T(1, 0) = (0, 1)$; $T(0, 1) = (-1, 0)$.
       The matrix is $R = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$.

8. **[Idempotency] Prove: If $P$ is a projection operator ($P^2=P$), its eigenvalues can only be 0 or 1.**
   ??? success "Solution"
       Let $Pv = \lambda v$, then $P^2 v = \lambda^2 v$.
       From $P^2=P$, $\lambda^2 v = \lambda v \implies (\lambda^2-\lambda)v = 0$.
       For a non-zero eigenvector $v$, $\lambda^2-\lambda=0$, so $\lambda=0$ or $1$.

9. **[Similarity] Do similar matrices have the same characteristic polynomial?**
   ??? success "Solution"
       Yes. $\det(\lambda I - P^{-1}AP) = \det(P^{-1}(\lambda I - A)P) = \det(P^{-1}) \det(\lambda I - A) \det(P) = \det(\lambda I - A)$.

10. **[Rank-Nullity App] What is the condition for $T: \mathbb{R}^n \to \mathbb{R}^m$ to be injective?**
    ??? success "Solution"
        The kernel must be zero, i.e., $\dim(\ker T) = 0$. By the rank-nullity theorem, this is equivalent to $\operatorname{rank}(T) = n$ (full column rank).

## Chapter Summary

Linear transformations unify operators and matrices:

1. **Representation**: A matrix is the numerical expression of a linear transformation in a specific coordinate system (basis).
2. **Internal Logic**: Similarity transformations reflect the description of the same operator from different viewpoints; its characteristic attributes remain invariant.
3. **Classification**: The dimensional distribution of kernel and image defines how much information is preserved or compressed by the transformation.
