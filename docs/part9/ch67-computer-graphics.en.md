# Chapter 67: Linear Algebra in Computer Graphics

<div class="context-flow" markdown>

**Prerequisites**: Linear Transformations (Ch05) · Homogeneous Coordinates · Geometric Algebra (Ch50) · Projections (Ch07)

**Chapter Outline**: From Geometric Entities to Pixel Mappings → Algebraic Construction of Homogeneous Coordinates → Unified Representation of Affine Transforms: Translation, Rotation, and Scaling → View Transformations and Perspective Projection → Depth Buffering and Z-transformations → Normal Vector Transformation via the Inverse-Transpose Matrix → Linear Control Point Representation of Splines → Applications: 3D Game Engine Architectures, Skeletal Animation in Movies, Ray Tracing, and Virtual Reality (VR)

**Extension**: Computer graphics is the "visual manifestation" of linear algebra; it maps rigid body motion and light propagation in the 3D world into a series of concatenated 4x4 matrix multiplications. It proves that realistic visual experiences are essentially matrix operations in homogeneous projective geometry—the mathematical soul of modern visual computing power.

</div>

In 3D games or movies, every frame on the screen is generated through millions of matrix multiplications. The core task of **Computer Graphics** is to project virtual 3D models onto a 2D screen. **Linear Algebra** provides the necessary tools: it uses homogeneous coordinates to unify translation and rotation into matrix operations and employs projection matrices to simulate the perspective effects of human vision. This chapter introduces the algebraic techniques serving as the bedrock of virtual reality.

---

## 67.1 Homogeneous Coordinates and Affine Transforms

!!! definition "Definition 67.1 (Homogeneous Coordinates)"
    To unify translation with linear operations, a 3D point $(x, y, z)$ is extended to a 4D vector $(x, y, z, 1)^T$.
    - **Translation Matrix $T(\mathbf{d})$**: $\begin{pmatrix} I_3 & \mathbf{d} \\ 0 & 1 \end{pmatrix}$.
    - **Rotation Matrix $R$**: $\begin{pmatrix} R_{3\times 3} & 0 \\ 0 & 1 \end{pmatrix}$.

---

## 67.2 Perspective Projection Matrices

!!! technique "Technique: Perspective Projection"
    Perspective projection uses the principle of similar triangles to make distant objects appear smaller. Its standard matrix form is:
    $$P = \begin{pmatrix} \frac{1}{\tan(\alpha/2) \cdot aspect} & 0 & 0 & 0 \\ 0 & \frac{1}{\tan(\alpha/2)} & 0 & 0 \\ 0 & 0 & -\frac{f+n}{f-n} & -\frac{2fn}{f-n} \\ 0 & 0 & -1 & 0 \end{pmatrix}$$
    The final "w-divide" step achieves the non-linear perspective effect.

---

## 67.3 Transformation of Normal Vectors

!!! theorem "Theorem 67.1 (Normal Transformation Rule)"
    If a surface undergoes a non-uniform scaling transformation $M$, its normal vector $\mathbf{n}$ must be transformed using the **Inverse-Transpose Matrix** $(M^{-1})^T$ to remain perpendicular to the tangent plane.

---

## Exercises

**1. [Basics] Write the $4 \times 4$ homogeneous rotation matrix for a 90-degree rotation around the $z$-axis.**

??? success "Solution"
    **Construction:**
    1. The 3D rotation part: $R_z(90^\circ) = \begin{pmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}$.
    2. Embed in a $4 \times 4$ matrix:
    **Conclusion**: $R = \begin{pmatrix} 0 & -1 & 0 & 0 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$.

**2. [Calculation] Translate the point $P = (1, 2, 3)^T$ by the vector $\mathbf{d} = (1, 1, 1)^T$ using a homogeneous matrix.**

??? success "Solution"
    **Matrix Multiplication:**
    $\begin{pmatrix} 1 & 0 & 0 & 1 \\ 0 & 1 & 0 & 1 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} 1 \\ 2 \\ 3 \\ 1 \end{pmatrix} = \begin{pmatrix} 1(1) + 1(1) \\ 1(2) + 1(1) \\ 1(3) + 1(1) \\ 1 \end{pmatrix} = \begin{pmatrix} 2 \\ 3 \\ 4 \\ 1 \end{pmatrix}$.
    **Conclusion**: The new coordinates are $(2, 3, 4)$.

**3. [Scaling] If a transformation matrix is $\operatorname{diag}(2, 2, 2, 1)$, how does the volume of the object change?**

??? success "Solution"
    **Analysis:**
    1. This is a uniform scaling factor of 2 along all three axes.
    2. The determinant $\det(M_{3\times 3}) = 2^3 = 8$.
    **Conclusion**: The volume of the object increases 8-fold.

**4. [Commutativity] In graphics, why is $TR \neq RT$?**

??? success "Solution"
    **Geometric Logic:**
    - $TR$: Rotate first, then translate in the world coordinate system.
    - $RT$: Translate first, then rotate around the world origin.
    Since matrix multiplication is non-commutative, the order determines whether the object rotates "in place" or "swings in an arc" around the origin.

**5. [Projection] What is the role of the last component $w$ in homogeneous coordinates?**

??? success "Solution"
    **Multiple Roles:**
    1. **Distinguish Vector vs Point**: $w=1$ for points, $w=0$ for vectors (translation has no effect on directions).
    2. **Enable Perspective**: After perspective transformation, $w$ becomes $-z$. By dividing all coordinates by $w$, the non-linear "near-large, far-small" mapping is realized.

**6. [Normals] Why can't normal vectors be multiplied directly by the transformation matrix $M$?**

??? success "Solution"
    **Reasoning:**
    Consider non-uniform scaling, e.g., a stretch by 2 along the $x$-axis. If a tangent vector is $(1, 1)$, the normal is $(-1, 1)$. After transformation, the tangent becomes $(2, 1)$. If the normal were multiplied by $M$, it would become $(-2, 1)$. Their dot product $(2)(-2) + (1)(1) = -3 \neq 0$; they are no longer perpendicular. Using $(M^{-1})^T$ yields the normal $(-0.5, 1)$, whose dot product with $(2, 1)$ is 0.

**7. [View Matrix] What is the "View Matrix" in graphics?**

??? success "Solution"
    It is the matrix that transforms objects from the "World Coordinate System" to the "Camera Coordinate System." It is essentially the **inverse** of the camera's own pose matrix. Linear algebra proves that "moving the camera forward" is mathematically identical to "moving the entire world backward."

**8. [Calculation] Convert the homogeneous coordinates $(4, 2, 8, 2)^T$ back to standard 3D coordinates.**

??? success "Solution"
    **Perform Perspective Divide:**
    $x = 4/2 = 2$
    $y = 2/2 = 1$
    $z = 8/2 = 4$
    **Conclusion**: The 3D coordinates are $(2, 1, 4)$.

**9. [Animation] Briefly explain the idea of Linear Blend Skinning (LBS).**

??? success "Solution"
    A vertex $v$ on an avatar's skin is influenced by adjacent bones (represented by matrices $M_1, M_2$). The final position is $v' = (\alpha M_1 + (1-\alpha) M_2) v$. This utilizes the **linear combination** of matrix space to allow skin to bend smoothly with bones.

**10. [Application] Why are GPUs designed as parallel matrix calculation architectures?**

??? success "Solution"
    **Reasoning:**
    Rendering involves performing the exact same $4 \times 4$ matrix multiplication on millions of independent vertices. This task lacks logical dependency, making it a classic Single-Instruction Multiple-Data (SIMD) problem. By executing these matrix operations across thousands of cores, GPUs reduce rendering time from minutes to milliseconds.

## Chapter Summary

Linear algebra serves as the "laws of physics" for virtual worlds:

1.  **Homogeneous Harmony**: Homogeneous coordinates unify complex Euclidean transformations into simple matrix multiplications, establishing the standard computation protocol for the graphics pipeline.
2.  **Dimensional Projection**: The perspective projection matrix proves that depth perception is essentially a non-linear projection operator, revealing the algebraic essence of optical imaging.
3.  **The Art of Computation**: From the inverse-transpose of normals to linear blend skinning, every realistic detail in graphics stems from precise control over the geometric properties of matrix space, supporting the visual revolution of modern digital media.
