# Chapter 47A: Matrix Calculus: Elementary

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Matrix Analysis (Ch14) · Multivariate Calculus Basics

**Chapter Outline**: From Scalar to Matrix Derivatives → Layout Conventions: Numerator vs. Denominator Layout → Gradients of Scalars with respect to Vectors → Gradients of Scalars with respect to Matrices → Key Formulas: Gradients of Linear Forms $a^T x$ and Quadratic Forms $x^T Ax$ → Derivatives of Trace Functions (Trace Tricks) → Matrix Form of the Chain Rule → Applications: Deriving Ordinary Least Squares (OLS), Gradients for Logistic Regression, and Basics of Linear Neural Networks

**Extension**: Matrix calculus is the symbolic engine of modern optimization algorithms; it simplifies tedious component-wise differentiation into concise matrix multiplications, serving as the essential mathematical path for understanding Backpropagation and any gradient-based machine learning model.

</div>

In multi-dimensional optimization problems, we frequently seek directions that maximize or minimize a scalar function (such as error, energy, or cost). **Matrix Calculus** provides a highly compact symbolic system for handling these problems. It not only significantly reduces the likelihood of manual derivation errors but also directly reveals the algebraic structure of the optimal solution. This chapter establishes standard differentiation rules and introduces "Trace Tricks" commonly used in statistics and machine learning.

---

## 47A.1 Layout Conventions and Basic Derivatives

!!! definition "Definition 47A.1 (Gradient Vector)"
    For a scalar function $f(\mathbf{x})$, its **gradient** with respect to an $n$-dimensional vector $\mathbf{x}$ is defined as:
    - **Denominator Layout** (Commonly used): $\frac{\partial f}{\partial \mathbf{x}}$ is an $n \times 1$ column vector.
    - **Numerator Layout**: $\frac{\partial f}{\partial \mathbf{x}}$ is a $1 \times n$ row vector.
    This encyclopedia defaults to the **Denominator Layout**.

!!! theorem "Theorem 47A.1 (Linear and Quadratic Gradients)"
    1.  $\frac{\partial (\mathbf{a}^T \mathbf{x})}{\partial \mathbf{x}} = \mathbf{a}$
    2.  $\frac{\partial (\mathbf{x}^T A \mathbf{x})}{\partial \mathbf{x}} = (A + A^T)\mathbf{x}$ (reduces to $2A\mathbf{x}$ if $A$ is symmetric).

---

## 47A.2 Trace Tricks

!!! technique "Technique: Derivatives of Trace Functions"
    For scalar functions expressed as traces $f(X) = \operatorname{tr}(AXB)$, the following formulas are foundational:
    1.  $\frac{\partial \operatorname{tr}(AX)}{\partial X} = A^T$
    2.  $\frac{\partial \operatorname{tr}(X^T A X)}{\partial X} = (A + A^T)X$
    Since any scalar $s$ satisfies $s = \operatorname{tr}(s)$, many complex quadratic derivatives can be solved easily by converting them into trace forms.

---

## 47A.3 Chain Rule and Hessians

!!! definition "Definition 47A.2 (Hessian Matrix)"
    The matrix of second-order partial derivatives of a scalar function $f$ with respect to a vector $\mathbf{x}$ is the **Hessian Matrix**:
    $$\mathbf{H} = \frac{\partial^2 f}{\partial \mathbf{x} \partial \mathbf{x}^T}$$
    **Significance**: The Hessian characterizes the curvature of the function; its positive definiteness determines whether a stationary point is a local minimum or maximum.

---

## Exercises

**1. [Basics] Calculate the gradient of $f(x) = \mathbf{x}^T \mathbf{x}$ with respect to $\mathbf{x}$.**

??? success "Solution"
    **Steps:**
    1. Express as a quadratic form: $f(x) = \mathbf{x}^T I \mathbf{x}$.
    2. Apply the quadratic gradient formula: $\frac{\partial f}{\partial \mathbf{x}} = (I + I^T)\mathbf{x} = 2I\mathbf{x}$.
    **Conclusion**: $\nabla f = 2\mathbf{x}$. This matches the scalar derivative of $x^2$ being $2x$.

**2. [Trace] Find $\frac{\partial \operatorname{tr}(A X^T)}{\partial X}$.**

??? success "Solution"
    **Steps:**
    1. Use trace property: $\operatorname{tr}(A X^T) = \operatorname{tr}(X A^T)$.
    2. Apply basic formula $\frac{\partial \operatorname{tr}(XB)}{\partial X} = B^T$.
    3. Set $B = A^T$.
    **Conclusion**: The derivative is $(A^T)^T = A$.

**3. [Linear Regression] Derive the gradient of the OLS loss function $J(\mathbf{w}) = \|\mathbf{y} - X\mathbf{w}\|^2$ with respect to $\mathbf{w}$.**

??? success "Solution"
    **Derivation:**
    1. Expand the norm: $J = (\mathbf{y}-X\mathbf{w})^T(\mathbf{y}-X\mathbf{w}) = \mathbf{y}^T\mathbf{y} - 2\mathbf{y}^TX\mathbf{w} + \mathbf{w}^TX^TX\mathbf{w}$.
    2. Differentiate with respect to $\mathbf{w}$:
       - Term 1 is constant: 0.
       - Term 2 is linear: $-2X^T\mathbf{y}$.
       - Term 3 is quadratic with $A=X^TX$: $2X^TX\mathbf{w}$.
    **Conclusion**: $\nabla_{\mathbf{w}} J = 2X^TX\mathbf{w} - 2X^T\mathbf{y}$. Setting this to zero yields the normal equation $X^TX\mathbf{w} = X^T\mathbf{y}$.

**4. [Hessian] Calculate the Hessian of $J(\mathbf{w})$ from the previous problem.**

??? success "Solution"
    **Steps:**
    1. Differentiate the first-order gradient $\mathbf{g} = 2X^TX\mathbf{w} - 2X^T\mathbf{y}$ again.
    2. The first term is linear in $\mathbf{w}$ with coefficient $2X^TX$.
    3. The second term is constant with respect to $\mathbf{w}$.
    **Conclusion**: $H = 2X^TX$. Since $X^TX$ is always positive semi-definite, this explains why the loss surface of least squares is a convex paraboloid.

**5. [Determinant] Given $\frac{\partial \ln \det X}{\partial X} = (X^{-1})^T$ for $X \succ 0$, find $\frac{\partial \det X}{\partial X}$.**

??? success "Solution"
    **Using Chain Rule:**
    1. Let $y = \det X$, then $\ln y = \ln \det X$.
    2. $\frac{\partial \ln y}{\partial X} = \frac{1}{y} \frac{\partial y}{\partial X}$.
    3. We know $\frac{\partial \ln y}{\partial X} = (X^{-1})^T$.
    4. Therefore, $\frac{\partial \det X}{\partial X} = (\det X) (X^{-1})^T$.
    **Note**: This is the transpose of the adjugate matrix $A^*$.

**6. [Verification] If $f = \mathbf{a}^T X \mathbf{b}$, find $\frac{\partial f}{\partial X}$.**

??? success "Solution"
    **Trace Trick:**
    1. Since $f$ is a scalar, $f = \operatorname{tr}(\mathbf{a}^T X \mathbf{b})$.
    2. Use cyclic property: $= \operatorname{tr}(\mathbf{b} \mathbf{a}^T X)$.
    3. Differentiate with respect to $X$ to get $(\mathbf{b} \mathbf{a}^T)^T = \mathbf{a} \mathbf{b}^T$.
    **Conclusion**: The gradient is a rank-1 matrix.

**7. [Frobenius] Compute $\frac{\partial \|X\|_F^2}{\partial X}$.**

??? success "Solution"
    **Calculation:**
    1. $\|X\|_F^2 = \operatorname{tr}(X^T X)$.
    2. Using the quadratic trace formula with $A=I$: $(I+I)X = 2X$.
    **Conclusion**: The derivative is $2X$.

**8. [Layout] In denominator layout, if $y = Ax$, what is the dimension of $\frac{\partial y}{\partial x}$?**

??? success "Solution"
    **Conclusion: $n \times m$ (specifically $A^T$).**
    In denominator layout, the gradient matrix of a vector with respect to a vector is the transpose of the Jacobian.

**9. [Non-linear] Find $\frac{\partial \exp(\mathbf{a}^T \mathbf{x})}{\partial \mathbf{x}}$.**

??? success "Solution"
    **Chain Rule:**
    1. Let $u = \mathbf{a}^T \mathbf{x}$.
    2. Then $f = e^u$.
    3. $\frac{\partial f}{\partial \mathbf{x}} = \frac{df}{du} \frac{\partial u}{\partial \mathbf{x}} = e^u \cdot \mathbf{a}$.
    **Conclusion**: $\exp(\mathbf{a}^T \mathbf{x}) \mathbf{a}$.

**10. [Application] Why is matrix calculus indispensable for weight updates in Deep Learning?**

??? success "Solution"
    **Reasoning:**
    Neural networks contain millions of parameters, usually arranged as weight matrices $W$. To minimize loss $L$, gradients $\nabla_W L$ must be computed. Matrix calculus allows us to write concise update rules like $\nabla_W L = \delta \cdot \mathbf{x}^T$ directly, rather than differentiating every single neuron connection individually, vastly improving implementation efficiency and speed.

## Chapter Summary

Matrix calculus transforms linear algebra into a symbolic machine for optimization:

1.  **Compression of Symbols**: Layout conventions and trace tricks compress thousands of partial derivatives into simple matrix products, enabling the symbolization of calculus in high-dimensional spaces.
2.  **Analytic Optimality**: Optimal solutions for models like linear regression are the algebraic roots of matrix gradient equations $\nabla f = 0$, revealing the necessary link between error surfaces and operator inverses.
3.  **Engine of Intelligence**: As the foundation for modern gradient descent, matrix calculus establishes the mathematical laws governing parameter flows in perceptrons, logistic regression, and deep neural networks.
