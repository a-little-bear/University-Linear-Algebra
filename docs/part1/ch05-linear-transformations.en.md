# Chapter 05: Linear Transformations

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch4) · Matrix Algebra (Ch2) · Function Theory

**Chapter Outline**: Definition of Linear Transformations → Kernel and Image (Range) → Rank-Nullity Theorem → Matrix Representation of a Transformation → Similarity and Coordinate Changes → Injective, Surjective, and Bijective Maps → Isomorphisms → Projections, Reflections, and Rotations → Composition of Transformations

**Extension**: Linear transformations are the "actions" of linear algebra; they represent how space is distorted, projected, or rotated while preserving the origin and straight lines.

</div>

Linear transformations are the mappings that preserve the algebraic structure of vector spaces. They are the functions that satisfy $T(u+v) = T(u) + T(v)$ and $T(ku) = k T(u)$. By fixing bases, every linear transformation between finite-dimensional spaces can be represented as a matrix. This chapter bridges the gap between abstract operators and concrete matrix calculations, culminating in the **Rank-Nullity Theorem**, which links the dimension of the input space to the dimensions of the kernel and range.

---

## 05.1 Definitions and Fundamental Theorem

!!! definition "Definition 05.1 (Linear Transformation)"
    A function $T: V \to W$ is a linear transformation if for all $u, v \in V$ and $c \in F$:
    1. $T(u + v) = T(u) + T(v)$
    2. $T(cu) = c T(u)$

!!! theorem "Theorem 05.1 (Rank-Nullity Theorem)"
    For a linear transformation $T: V \to W$:
    $$\dim(\ker T) + \dim(\operatorname{Im} T) = \dim V$$
    The dimension of the kernel plus the dimension of the image equals the dimension of the domain.

---

## Exercises

1. **[Fundamentals] Is the map $T(x, y) = (x+1, y)$ a linear transformation?**
   ??? success "Solution"
       No. $T(0, 0) = (1, 0) \neq (0, 0)$. A linear transformation must map the zero vector to the zero vector.

2. **[Kernel] Find the kernel of the projection map $P(x, y, z) = (x, y, 0)$.**
   ??? success "Solution"
       The kernel consists of all vectors $(x, y, z)$ such that $(x, y, 0) = (0, 0, 0)$. This implies $x=0$ and $y=0$. Thus $\ker P = \{ (0, 0, z) : z \in \mathbb{R} \}$, which is the $z$-axis.

3. **[Image] Describe the image of the differentiation operator $D: P_3 \to P_2$ where $D(p) = p'$.**
   ??? success "Solution"
       The image consists of all derivatives of polynomials of degree at most 3. Since every polynomial of degree 2 is the derivative of some polynomial of degree 3, $\operatorname{Im} D = P_2$.

4. **[Matrix Representation] Find the matrix of $T(x, y) = (2x-y, x+3y)$ relative to the standard basis.**
   ??? success "Solution"
       $T(1, 0) = (2, 1)$ and $T(0, 1) = (-1, 3)$. The matrix is $\begin{pmatrix} 2 & -1 \\ 1 & 3 \end{pmatrix}$.

5. **[Similarity] When are two matrices $A$ and $B$ called similar?**
   ??? success "Solution"
       $A$ and $B$ are similar if there exists an invertible matrix $P$ such that $B = P^{-1} A P$. Similar matrices represent the same linear transformation in different bases.

6. **[Isomorphism] Define a linear isomorphism.**
   ??? success "Solution"
       A linear isomorphism is a bijective linear transformation. If such a map exists between $V$ and $W$, the spaces are structurally identical as vector spaces.

7. **[Rank-Nullity] If $T: \mathbb{R}^5 \to \mathbb{R}^3$ is surjective, what is $\dim \ker T$?**
   ??? success "Solution"
       Surjective implies $\dim(\operatorname{Im} T) = 3$. By Rank-Nullity: $\dim \ker T + 3 = 5$, so $\dim \ker T = 2$.

8. **[Injective] Show that $T$ is injective iff $\ker T = \{0\}$.**
   ??? success "Solution"
       If $\ker T = \{0\}$, then $T(u) = T(v) \implies T(u-v) = 0 \implies u-v = 0 \implies u = v$. Conversely, if $T$ is injective, $T(v) = 0 = T(0)$ implies $v=0$.

9. **[Projection] Prove that for any projection operator $P^2 = P$, $V = \ker P \oplus \operatorname{Im} P$.**
   ??? success "Solution"
       For any $v$, $v = (v - Pv) + Pv$. Note $P(v - Pv) = Pv - P^2v = Pv - Pv = 0$, so $(v - Pv) \in \ker P$. Also $Pv \in \operatorname{Im} P$. The intersection is $\{0\}$ since if $Pv = 0$ and $v = Pu$, then $v = Pu = P^2u = P(Pu) = Pv = 0$.

10. **[Trace] Why is $\operatorname{tr}(A) = \operatorname{tr}(P^{-1}AP)$?**
    ??? success "Solution"
        Using the property $\operatorname{tr}(XY) = \operatorname{tr}(YX)$, we have $\operatorname{tr}(P^{-1}(AP)) = \operatorname{tr}((AP)P^{-1}) = \operatorname{tr}(A)$. Thus the trace is a well-defined invariant of the linear transformation.

## Chapter Summary

This chapter explores the functional perspective of linear algebra:

1. **Structural Preservation**: Defined linear transformations as mappings that respect vector addition and scaling.
2. **Space Partitioning**: Developed the kernel and image as the fundamental descriptors of an operator's behavior.
3. **Dimensional Conservation**: Established the Rank-Nullity Theorem as the governing law for the mapping of dimensions.
4. **Coordinate Independence**: Linked transformations to matrices, identifying similarity as the equivalence relation for basis changes.
