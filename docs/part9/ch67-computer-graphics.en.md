# Chapter 67: Linear Algebra in Computer Graphics

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Linear Transformations (Ch05) · Quaternion Matrices (Ch51) · Orthogonality (Ch07)

**Chapter Outline**: Introduction of Homogeneous Coordinates → Unified Matrix Representation of Translation, Rotation, and Scaling → 2D and 3D Transformation Matrices → Projection Transforms: Orthographic and Perspective → View and Camera Models → Barycentric Coordinates & Triangle Rasterization → Linear Intersection Equations in Ray Tracing → Shading Models: Dot Products and Lambertian Reflectance → Animation Interpolation: LERP and SLERP

**Extension**: Computer graphics is essentially "real-time visualization" of linear algebra; it projects the geometry of a 3D world onto a 2D screen using 4x4 matrices, forming the mathematical heart of modern game engines (Unreal, Unity) and visual effects.

</div>

Every 3D character or scene you see on a computer screen is the result of millions of matrix multiplications. To handle translation (non-linear) and rotation (linear) within a unified framework, graphics utilizes **Homogeneous Coordinates**. This chapter demonstrates how to build the complete geometric path from a virtual world to a pixel plane using 4x4 matrices.

---

## 67.1 Homogeneous Coordinates and Basic Transforms

!!! definition "Definition 67.1 (Homogeneous Coordinates)"
    To represent translation as matrix multiplication, we augment a 3D vector $(x, y, z)$ into a 4D vector $(x, y, z, w)$, typically setting $w=1$.
    - **Translation Matrix**: $T(\mathbf{t}) = \begin{pmatrix} 1 & 0 & 0 & t_x \\ 0 & 1 & 0 & t_y \\ 0 & 0 & 1 & t_z \\ 0 & 0 & 0 & 1 \end{pmatrix}$.

!!! technique "Composition of Transformations"
    Using matrix multiplication, we can combine scaling ($S$), rotation ($R$), and translation ($T$) into a single transformation matrix $M = T \cdot R \cdot S$. Note that the order of multiplication determines the relative coordinate frame of each transform.

---

## 67.2 Projection Transformations

!!! definition "Definition 67.2 (Perspective Projection Matrix)"
    Perspective projection simulates the "smaller when farther" characteristic of human vision. Given frustum parameters, its standard matrix form is:
    $$P = \begin{pmatrix} \frac{1}{\tan(\alpha/2)} & 0 & 0 & 0 \\ 0 & \frac{1}{\tan(\alpha/2)} & 0 & 0 \\ 0 & 0 & \frac{f+n}{f-n} & \frac{-2fn}{f-n} \\ 0 & 0 & 1 & 0 \end{pmatrix}$$
    The last row copies the $z$ coordinate into the $w$ component to enable perspective division ($x/w, y/w$).

---

## 67.3 Barycentric Coordinates and Rasterization

!!! technique "Barycentric Coordinates"
    Any point $P$ inside a triangle can be expressed as a linear combination of its three vertices $A, B, C$:
    $$P = \alpha A + \beta B + \gamma C, \quad \alpha + \beta + \gamma = 1$$
    **Application**: This is the core algorithm for attribute interpolation (color, texture coordinates) and for testing whether a point lies inside a triangle.

---

## 67.4 Shading and Lighting

!!! theorem "Theorem 67.1 (Lambertian Reflectance)"
    The intensity $I$ of light at a point on a surface is proportional to the dot product of the surface normal $\mathbf{N}$ and the light direction $\mathbf{L}$:
    $$I = I_0 \max(0, \mathbf{N} \cdot \mathbf{L})$$
    This shows how the **Inner Product** in linear algebra determines shadows and the sense of 3D volume.

---

## Exercises

1.  **[Basics] Write the 2D homogeneous transformation matrix for rotation by angle $\theta$ around the origin.**
    ??? success "Solution"
        $\begin{pmatrix} \cos\theta & -\sin\theta & 0 \\ \sin\theta & \cos\theta & 0 \\ 0 & 0 & 1 \end{pmatrix}$.

2.  **[Translation] Find the homogeneous coordinates of the point $(1, 2, 3)$ translated by $(5, 0, -1)$.**
    ??? success "Solution"
        $(1+5, 2+0, 3-1, 1) = (6, 2, 2, 1)$.

3.  **[Projection] Why is the last row of a perspective projection matrix usually $(0, 0, 1, 0)$?**
    ??? success "Solution"
        To store the $z$-value in the homogeneous $w$-component. Hardware then performs the normalization $x/w, y/w, z/w$, creating the perspective effect.

4.  **[Barycentric] If a point has barycentric coordinates $(0.5, 0.5, 0)$, where is it located?**
    ??? success "Solution"
        It is at the midpoint of the edge connecting the first two vertices ($AB$).

5.  **[Quaternions] Why are quaternions used for 3D rotation in games instead of Euler angles?**
    ??? success "Solution"
        To avoid **Gimbal Lock**, and because spherical linear interpolation (SLERP) is smooth and more computationally efficient.

6.  **[Ray Tracing] Calculate the intersection $t$ of a ray $\mathbf{r}(t) = \mathbf{o} + t\mathbf{d}$ and a plane $\mathbf{n} \cdot \mathbf{p} + d = 0$.**
    ??? success "Solution"
        Substitute: $\mathbf{n} \cdot (\mathbf{o} + t\mathbf{d}) + d = 0 \implies t = \frac{-(d + \mathbf{n} \cdot \mathbf{o})}{\mathbf{n} \cdot \mathbf{d}}$.

7.  **[Order] Prove that "rotate then translate" is not equivalent to "translate then rotate."**
    ??? success "Solution"
        $T \cdot R \neq R \cdot T$. Matrix multiplication is non-commutative; translating and then rotating will cause the translation vector itself to rotate.

8.  **[Normals] Why must we transform normals $\mathbf{N}$ using $(M^{-1})^T$ instead of the model matrix $M$?**
    ??? success "Solution"
        To preserve the orthogonality between the normal and the tangent plane. Under non-uniform scaling, $M$ would distort the normal's direction.

9.  **[SLERP] What is Spherical Linear Interpolation?**
    ??? success "Solution"
        Constant angular velocity interpolation between two unit quaternions (points on a 4D sphere) along a great-circle path.

****

??? success "Solution"
    

## Chapter Summary

Computer graphics is the most spectacular stage for linear algebra:

1.  **Elevation of Dimension**: Through homogeneous coordinates, graphics unifies non-linear rigid body motion into linear matrix arithmetic, establishing the algebraic basis for pipelined rendering.
2.  **Mystery of Projection**: The perspective matrix demonstrates how to utilize the rank-deficiency tendencies of linear algebra (via $w$-division) to simulate the collapse of high-dimensional space into low-dimensional vision, constructing the sense of depth in digital worlds.
3.  **Logic of Geometry**: From barycentric interpolation to dot-product lighting, graphics proves that all real-world sensory experiences (color, shape, light) can be reduced to elegant vector operations at the underlying level.
