# Chapter 67: Applications of Linear Algebra in Computer Graphics

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · Linear Transformations (Ch5) · Orthogonal Matrices (Ch7) · Eigenvalues (Ch6)

**Chapter Outline**: Homogeneous Coordinates → Affine Transformations → Rotation Representations (Matrix/Quaternion/Euler/SLERP) → Model-View-Projection Matrix → Perspective and Orthographic Projections → Viewport Transformation → Color Space Transformations → Spline Surfaces

**Extension**: Computer graphics is a high-ground for applied linear algebra; GPU architectures are essentially large-scale Single Instruction Multiple Data (SIMD) matrix-multiplication accelerators.

</div>

Computer graphics achieves digital reconstruction of 3D scenes through linear space transformations. By introducing homogeneous coordinates, non-linear affine transformations are embedded into high-dimensional linear mappings, allowing complex geometric pipelines to be handled through unified matrix multiplication chains.

---

## 67.1 Homogeneous Transformations and Geometric Pipelines

!!! definition "Definition 67.1 (Homogeneous Coordinates and Affine Embedding)"
    A point $\mathbf{p} \in \mathbb{R}^3$ in 3D Euclidean space is mapped to a quadruple $(x, y, z, 1)^T$ in projective space. This allows translation transformations $T$ to be represented as elements of a subgroup of $GL(4, \mathbb{R})$, alongside rotations and scaling.

!!! theorem "Theorem 67.7 (Adjoint Normal Transformation)"
    If vertices are transformed by a non-singular matrix $M$, then their surface normals $\mathbf{n}$ must follow the transformation rule of $(M^{-1})^T$ to maintain orthogonality between the normal and the tangent plane.

---

## Exercises

1. **[Matrix Unity] Explain why homogeneous coordinates are the foundation for achieving standardized processing in GPU hardware pipelines.**
   ??? success "Solution"
       Homogeneous coordinates convert translation (addition) into matrix multiplication. This allows rotations, scaling, translations, and projections to be concatenated into a single composite matrix $M_{total}$. GPU hardware only needs to implement efficient matrix-vector multipliers to process all geometric operations.

2. **[Homogeneous Division] Derive the logic for generating the $w$ component in a perspective projection matrix and its role in Homogeneous Division.**
   ??? success "Solution"
       The perspective projection matrix encodes $z$ information into the $w$ component (usually $w = -z$). After clipping, all coordinates are divided by $w$, achieving the mapping of $x/z, y/z$. This is geometrically equivalent to similar triangle ratios, producing the visual perspective effect.

3. **[Rotational Singularity] Analyze the mathematical essence of "Gimbal Lock" in Euler angle representation when $	heta = \pm 90^\circ$.**
   ??? success "Solution"
       When the pitch angle is $\pm 90^\circ$, the axes for roll and yaw rotations become collinear. At this point, the Jacobian of the rotation matrix becomes rank-deficient, and the system loses a degree of freedom, making it impossible to describe small rotations in certain directions in 3D space.

4. **[Quaternions] Prove that a unit quaternion $\mathbf{q}$ and its negative $-\mathbf{q}$ describe the same rotation operator.**
   ??? success "Solution"
       The rotation transformation is given by $\mathbf{p}' = \mathbf{q} \mathbf{p} \mathbf{q}^{-1}$. Substituting $-\mathbf{q}$ yields $(-\mathbf{q}) \mathbf{p} (-\mathbf{q})^{-1} = (-1)^2 \mathbf{q} \mathbf{p} \mathbf{q}^{-1} = \mathbf{p}'$. This demonstrates that $SU(2)$ is a double cover of $SO(3)$.

5. **[Interpolation Properties] Prove that the SLERP (Spherical Linear Interpolation) formula has constant angular velocity on the unit sphere.**
   ??? success "Solution"
       The SLERP formula is given by $\frac{\sin((1-t)	heta)}{\sin	heta} q_0 + \frac{\sin(t	heta)}{\sin	heta} q_1$. By calculating the norm of its first derivative with respect to $t$, one can prove that its arc length evolution rate on $S^3$ is constant, i.e., $\frac{d\alpha}{dt} = 	ext{const}$.

6. **[Normal Correction] Calculate the effect of the non-uniform scaling matrix $S = \operatorname{diag}(2, 1, 1)$ on the normal vector $(1, 1, 0)^T$.**
   ??? success "Solution"
       The vertex matrix is $M=S$. The correct normal transformation matrix is $(M^{-1})^T = \operatorname{diag}(0.5, 1, 1)$. The transformed normal vector is $(0.5, 1, 0)^T$. To verify: consider a tangent vector $(1, -1, 0)^T$ (orthogonal to the original normal). The transformed tangent is $(2, -1, 0)^T$. The inner product is $0.5 \cdot 2 + 1 \cdot (-1) = 0$, confirming the normal remains orthogonal.

7. **[Depth Nonlinearity] Analyze the reciprocal relationship between the depth value $z_{ndc}$ after perspective transformation and the original depth $z_{eye}$, and its impact on Z-buffer precision distribution.**
   ??? success "Solution"
       $z_{ndc} = A + B/z_{eye}$. This hyperbolic relationship results in an extremely high sampling density near the near plane and very sparse sampling near the far plane. This favors precise occlusion determination for close objects but leads to depth conflicts (Z-fighting) in distant views.

8. **[Barycentric Coordinates] Prove that the perspective-correct interpolation formula $a = \frac{\sum a_i \lambda_i / w_i}{\sum \lambda_i / w_i}$ restores linearity of attributes in 3D space.**
   ??? success "Solution"
       This stems from the fact that the projective transformation from 3D world space to screen space is not affine. Only through $1/w$ weighting can the non-linear screen coordinate system be re-mapped back to the linearly varying clip space parameters.

9. **[Reflection and Involution] Prove that the Householder reflection matrix $H = I - 2\mathbf{n}\mathbf{n}^T$ is symmetric and orthogonal, and find its determinant.**
   ??? success "Solution"
       $H^T = H$ is obvious. $H^T H = (I-2\mathbf{n}\mathbf{n}^T)(I-2\mathbf{n}\mathbf{n}^T) = I - 4\mathbf{n}\mathbf{n}^T + 4\mathbf{n}(\mathbf{n}^T\mathbf{n})\mathbf{n}^T = I$. The determinant is calculated via $\det(I + uv^T) = 1 + v^T u$ as $1 + (-2)\mathbf{n}^T\mathbf{n} = -1$.

10. **[Color Transformation] Explain the basis for determining the luminance component $Y$ coefficients in the linear transformation matrix from RGB space to YCbCr space.**
    ??? success "Solution"
        $Y = 0.299R + 0.587G + 0.114B$. The coefficients reflect the weighted sensitivity of the human visual system (CIE observer model) to different wavelengths of light intensity (strongest for green, weakest for blue).

## Chapter Summary

This chapter systematically discusses the foundational role of linear algebra in computer graphics:

1. **Unified Representation**: Homogeneous coordinates unify geometric transformations as $GL(4, \mathbb{R})$ matrix operations, establishing the computational standard for graphics pipelines.
2. **Rotation Modeling**: Analyzed the algebraic properties of various rotation representations, proving the superiority of quaternions in continuous motion interpolation.
3. **Spatial Projection**: Revealed the non-linear nature of perspective transformation and its linear implementation under projective geometry.
4. **Attribute Geometry**: Established criteria for normal vector transformation and perspective-correct interpolation, ensuring physical and geometric consistency in rendering.
