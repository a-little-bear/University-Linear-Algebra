# Chapter 46B: Matrix Means

<div class="context-flow" markdown>

**Prerequisites**: Operator Monotone Functions (Ch46A) · Positive Definite Matrices (Ch16) · Matrix Analysis (Ch14)

**Chapter Outline**: Operatorization of the Mean Concept → Kubo-Ando Axioms (Monotonicity, Continuity, Transformer Invariance) → Core Mapping: One-to-one Correspondence between Means and Operator Monotone Functions → Classic Operator Means: Arithmetic ($A \nabla B$), Geometric ($A \# B$), and Harmonic ($A ! B$) → Operator Mean Inequalities → Explicit Formulas for Matrix Geometric Means → Applications: Quantum Metric Geometry (Fisher Information Metric), Parallel Impedance in Electrical Engineering, and Covariance Averaging in Signal Processing

**Extension**: Matrix mean theory is the "harmonizing tool" for operator non-commutativity; it constructs an axiomatic system via operator monotone functions, proving that even in non-commutative environments, we can define "intermediate" operators with excellent algebraic properties—the mathematical core of modern quantum metrology.

</div>

In scalar algebra, $(a+b)/2$ and $\sqrt{ab}$ are elementary concepts. However, when $A$ and $B$ are non-commuting matrices, a direct definition like $\sqrt{AB}$ results in something that is generally not even Hermitian. The theory of **Matrix Means**, pioneered by **Kubo and Ando**, provides a rigorous "averaging" scheme for positive definite matrices. This chapter introduces how to generate physically meaningful matrix means using operator monotone functions.

---

## 46B.1 Kubo-Ando Axiomatic Definition

!!! definition "Definition 46B.1 (Operator Mean $\sigma$)"
    A binary operation $(A, B) \mapsto A \sigma B$ on positive operators is an operator mean if it satisfies:
    1.  **Monotonicity**: If $A \preceq A'$ and $B \preceq B'$, then $A \sigma B \preceq A' \sigma B'$.
    2.  **Transformer Invariance**: $C(A \sigma B)C^* = (CAC^*) \sigma (CBC^*)$.
    3.  **Continuity**: It is continuous under monotonic sequences.
    4.  **Normalization**: $I \sigma I = I$.

---

## 46B.2 Core Correspondence

!!! theorem "Theorem 46B.1 (Means vs. Functions)"
    Every operator mean $\sigma$ corresponds uniquely to a positive operator monotone function $f$ on $(0, \infty)$:
    $$A \sigma B = A^{1/2} f(A^{-1/2} B A^{-1/2}) A^{1/2}$$
    - **Arithmetic Mean**: $f(t) = (1+t)/2 \Rightarrow A \nabla B = (A+B)/2$.
    - **Geometric Mean**: $f(t) = \sqrt{t} \Rightarrow A \# B = A^{1/2} (A^{-1/2} B A^{-1/2})^{1/2} A^{1/2}$.
    - **Harmonic Mean**: $f(t) = 2t/(1+t) \Rightarrow A ! B = 2(A^{-1} + B^{-1})^{-1}$.

---

## 46B.3 Operator Mean Inequalities

!!! theorem "Theorem 46B.2 (Operator A-G-H Inequality)"
    For any positive definite matrices $A, B$:
    $$A ! B \preceq A \# B \preceq A \nabla B$$
    **Significance**: This chain of inequalities preserves the order of scalar means at the operator level, establishing a hierarchy of "intermediate values" for matrices.

---

## Exercises

**1. [Basics] Calculate the geometric mean of $A = \operatorname{diag}(1, 1)$ and $B = \operatorname{diag}(4, 9)$.**

??? success "Solution"
    **Since $A$ and $B$ commute (diagonal):**
    1. The geometric mean reduces to the entry-wise scalar geometric mean.
    2. Component 1: $\sqrt{1 \cdot 4} = 2$.
    3. Component 2: $\sqrt{1 \cdot 9} = 3$.
    **Conclusion**: $A \# B = \operatorname{diag}(2, 3)$.

**2. [Formula] Verify $\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} \# \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. The product of these matrices is the zero matrix.
    2. Although the formula involves $A^{-1/2}$, for semi-definite cases the mean is defined as a limit.
    3. Since the subspaces spanned by the two operators are orthogonal, their "overlap" is zero.
    **Conclusion**: $A \# B = O$.

**3. [Property] Prove the symmetry of the geometric mean: $A \# B = B \# A$.**

??? success "Solution"
    **Proof:**
    1. Formula: $A \# B = A^{1/2}(A^{-1/2} B A^{-1/2})^{1/2} A^{1/2}$.
    2. Using the identity $X(X^{-1}Y)^{1/2} = (YX^{-1})^{1/2}X$ (where $X=A^{1/2}, Y=B^{1/2}$), the expression can be rearranged into the symmetric form.
    **Conclusion**: The geometric mean is symmetric with respect to its arguments.

**4. [Harmonic] If $A$ and $B$ are resistance matrices, what physical configuration does $A ! B$ represent?**

??? success "Solution"
    **Physical Logic:**
    1. The total resistance of two resistors $R_1, R_2$ in parallel is $R_{total} = (R_1^{-1} + R_2^{-1})^{-1}$.
    2. The operator harmonic mean is $A ! B = 2(A^{-1} + B^{-1})^{-1}$.
    **Conclusion**: Excluding the factor of 2, this is the mathematical model for **parallel impedance** in multi-port networks.

**5. [Monotonicity] If $A \preceq A'$, prove $A \nabla B \preceq A' \nabla B$.**

??? success "Solution"
    **Proof:**
    1. $A \nabla B = (A+B)/2$.
    2. Difference: $(A'+B)/2 - (A+B)/2 = (A'-A)/2$.
    3. Since $A \preceq A'$, $A'-A \succeq 0$.
    4. Thus the difference is positive semi-definite. The property holds.

**6. [Calculation] Find the difference between the arithmetic and geometric means of $A = 2I$ and $B = 8I$.**

??? success "Solution"
    **Calculation:**
    1. $A \nabla B = (2+8)/2 \cdot I = 5I$.
    2. $A \# B = \sqrt{2 \cdot 8} \cdot I = 4I$.
    3. $A \nabla B - A \# B = I \succeq 0$. Consistently follows the A-G inequality.

**7. [Uniqueness] Why can't $A \# B$ be simply defined as $\sqrt{AB}$?**

??? success "Solution"
    **Reasoning:**
    1. If $A$ and $B$ do not commute, $AB$ is generally **not Hermitian** (symmetric).
    2. An operator mean must be defined within the space of Hermitian operators to maintain physical meaning and partial order properties.
    3. The "sandwich structure" $A^{1/2}(\cdot)A^{1/2}$ ensures the result is always symmetric and positive definite via a congruence transform.

**8. [Application] Briefly state the role of the matrix geometric mean in Diffusion Tensor Imaging (DTI).**

??? success "Solution"
    In DTI, each pixel represents a positive definite matrix (tensor). Using arithmetic means for smoothing can lead to "swelling" effects (increased volume). The geometric mean preserves determinant properties (volume) while providing a more natural interpolation along the manifold of positive definite matrices.

**9. [Commutative] Prove: If $A$ and $B$ commute, $A \# B = \sqrt{AB}$.**

??? success "Solution"
    **Proof:**
    1. If they commute, they share an eigenbasis.
    2. $A^{1/2} (A^{-1/2} B A^{-1/2})^{1/2} A^{1/2} = A^{1/2} (A^{-1} B)^{1/2} A^{1/2}$.
    3. $= A^{1/2} A^{-1/2} B^{1/2} A^{1/2} = B^{1/2} A^{1/2} = (BA)^{1/2} = (AB)^{1/2}$.

**10. [Maximal] Among all operator means, which is the largest?**

??? success "Solution"
    **Conclusion: The Arithmetic Mean $A \nabla B$.**
    **Reasoning**: Any operator monotone function $f$ with $f(1)=1$ satisfies $f(t) \le (1+t)/2$ on $(0, \infty)$. This ensures the arithmetic mean is the upper bound for any axiomatic mean system.

## Chapter Summary

Matrix mean theory provides the "middle way algebra" for operator interactions:

1.  **Axiomatic Balance**: The Kubo-Ando axioms define the boundaries of legitimacy for means, transforming intuitive feelings of "averaging" into rigorous traits like monotonicity and transformer invariance.
2.  **Non-commutative Bridge**: Via the "sandwich structure," concepts like the geometric mean successfully navigate the obstacles of non-commutativity, introducing natural geodesic paths to the space of positive operators.
3.  **Extending Inequalities**: The operator A-G-H inequalities are not just generalizations of scalar results but describe deep energy laws governing quantum information flow and physical impedance superposition.
