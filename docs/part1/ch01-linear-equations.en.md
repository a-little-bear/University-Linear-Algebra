# Chapter 1  Systems of Linear Equations

<div class="context-flow" markdown>

**Prerequisites**: None (starting point of linear algebra) · **Chapter arc**: Linear equations → augmented matrices → elementary row operations → Gaussian elimination → REF/RREF → existence and uniqueness of solutions → homogeneous systems → solution structure
**One-line essence**: The prototype of every linear algebra problem is $A\mathbf{x}=\mathbf{b}$ — this chapter establishes the complete algorithmic framework for solving it and the structural theory of solution sets

</div>

A system of linear equations is one of the most fundamental objects of study in linear algebra. From solving systems of two linear equations in elementary mathematics, we have already been dealing with problems of solving linear systems. In this chapter, we will systematically establish the theoretical framework for systems of linear equations: introduce augmented matrices and elementary row operations, develop Gaussian elimination as a powerful algorithmic tool, and deeply analyze the existence, uniqueness, and structure of solutions. These topics form the foundation for all subsequent chapters.

---

## 1.1 Linear Equations and Systems of Linear Equations

<div class="context-flow" markdown>

**Core question**: What is "linear"? → unknowns appear to the first power, no products → geometrically corresponds to the intersection of lines/planes/hyperplanes

</div>

### Basic Definitions

!!! definition "Definition 1.1 (Linear equation)"
    A **linear equation** in $n$ unknowns $x_1, x_2, \ldots, x_n$ is an equation of the form

    $$
    a_1 x_1 + a_2 x_2 + \cdots + a_n x_n = b
    $$

    where $a_1, a_2, \ldots, a_n$ and $b$ are known real (or complex) constants, called the **coefficients** and the **constant term** of the equation, respectively.

The key characteristic of a linear equation is that each unknown appears only to the first power, with no product terms between unknowns. For example, $3x_1 - 2x_2 + x_3 = 7$ is a linear equation, while $x_1 x_2 + x_3 = 1$ or $\sin(x_1) + x_2 = 0$ are not.

!!! definition "Definition 1.2 (System of linear equations)"
    A collection of $m$ linear equations in $n$ unknowns is called an $m \times n$ **system of linear equations**, with the general form:

    $$
    \begin{cases}
    a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n = b_1 \\
    a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n = b_2 \\
    \quad \vdots \\
    a_{m1}x_1 + a_{m2}x_2 + \cdots + a_{mn}x_n = b_m
    \end{cases}
    $$

    where $a_{ij}$ ($1 \le i \le m, 1 \le j \le n$) are the coefficients and $b_i$ are the constant terms.

!!! definition "Definition 1.3 (Solution and solution set)"
    A set of numbers $(s_1, s_2, \ldots, s_n)$ that simultaneously satisfies all equations in a system of linear equations is called a **solution** of the system. The set of all solutions is called the **solution set**. If two systems of linear equations have the same solution set, they are said to be **equivalent**.

### Geometric Interpretation

**Two-dimensional case (two unknowns)**: Each linear equation $a_1 x_1 + a_2 x_2 = b$ represents a line in the $\mathbb{R}^2$ plane. The solution set of two equations corresponds to the intersection of two lines, with three possible cases:

1. **Unique solution**: The two lines intersect at a single point (the lines have different slopes).
2. **Infinitely many solutions**: The two lines coincide (the equations are proportional).
3. **No solution**: The two lines are parallel but not coincident.

**Three-dimensional case (three unknowns)**: Each linear equation $a_1 x_1 + a_2 x_2 + a_3 x_3 = b$ represents a plane in $\mathbb{R}^3$. The solution set of three equations is the common intersection of three planes, which may be a point, a line, a plane, or the empty set.

!!! example "Example 1.1"
    Consider the system

    $$
    \begin{cases} x_1 + x_2 = 3 \\ 2x_1 - x_2 = 0 \end{cases}
    $$

    The first equation represents a line passing through $(3,0)$ and $(0,3)$ in the plane, and the second equation represents a line through the origin with slope $2$. The two lines intersect at the point $(1, 2)$, which is the unique solution $x_1 = 1, x_2 = 2$.

!!! example "Example 1.2"
    Consider the three-dimensional system

    $$
    \begin{cases} x_1 + x_2 + x_3 = 1 \\ x_1 - x_2 + x_3 = 3 \end{cases}
    $$

    Each equation represents a plane in $\mathbb{R}^3$. Two non-parallel planes intersect in a line, so this system has infinitely many solutions. Subtracting the two equations gives $2x_2 = -2$, i.e., $x_2 = -1$. Substituting into the first equation gives $x_1 + x_3 = 2$. Setting $x_3 = t$ (free parameter), the general solution is $x_1 = 2 - t,\; x_2 = -1,\; x_3 = t$, $t \in \mathbb{R}$.

---

## 1.2 Augmented Matrices and Elementary Row Operations

<div class="context-flow" markdown>

**Motivation**: The unknown symbols in systems are redundant → compress the information into an **augmented matrix** → use **elementary row operations** to systematically manipulate it (→ Chapter 2 abstracts row operations as elementary matrices)

</div>

When solving systems of linear equations, the symbols $x_1, x_2, \ldots, x_n$ for the unknowns serve only as placeholders; the coefficients and constant terms are what truly determine the solutions. Therefore, we can compress all the information of a system into a matrix.

!!! definition "Definition 1.4 (Coefficient matrix and augmented matrix)"
    For the system of linear equations

    $$
    \begin{cases}
    a_{11}x_1 + a_{12}x_2 + \cdots + a_{1n}x_n = b_1 \\
    a_{21}x_1 + a_{22}x_2 + \cdots + a_{2n}x_n = b_2 \\
    \quad \vdots \\
    a_{m1}x_1 + a_{m2}x_2 + \cdots + a_{mn}x_n = b_m
    \end{cases}
    $$

    its **coefficient matrix** is

    $$
    A = \begin{pmatrix}
    a_{11} & a_{12} & \cdots & a_{1n} \\
    a_{21} & a_{22} & \cdots & a_{2n} \\
    \vdots & \vdots & \ddots & \vdots \\
    a_{m1} & a_{m2} & \cdots & a_{mn}
    \end{pmatrix}
    $$

    its **augmented matrix** is

    $$
    [A \mid \mathbf{b}] = \left(\begin{array}{cccc|c}
    a_{11} & a_{12} & \cdots & a_{1n} & b_1 \\
    a_{21} & a_{22} & \cdots & a_{2n} & b_2 \\
    \vdots & \vdots & \ddots & \vdots & \vdots \\
    a_{m1} & a_{m2} & \cdots & a_{mn} & b_m
    \end{array}\right)
    $$

The system can also be written compactly as the matrix equation $A\mathbf{x} = \mathbf{b}$, where $\mathbf{x} = (x_1, x_2, \ldots, x_n)^T$, $\mathbf{b} = (b_1, b_2, \ldots, b_m)^T$.

!!! definition "Definition 1.5 (Elementary row operations)"
    The following three operations on a matrix are called **elementary row operations**:

    1. **Row interchange**: Swap rows $i$ and $j$ of the matrix, denoted $R_i \leftrightarrow R_j$.
    2. **Row scaling**: Multiply row $i$ by a nonzero constant $c$, denoted $cR_i \to R_i$.
    3. **Row replacement**: Add $c$ times row $j$ to row $i$, denoted $R_i + cR_j \to R_i$.

!!! theorem "Theorem 1.1 (Row operations preserve solution sets)"
    Applying finitely many elementary row operations to the augmented matrix of a system of linear equations produces a new system equivalent to the original, i.e., they have the same solution set.

??? proof "Proof"
    Each elementary row operation is reversible:

    - The inverse of row interchange $R_i \leftrightarrow R_j$ is again $R_i \leftrightarrow R_j$.
    - The inverse of row scaling $cR_i \to R_i$ ($c \neq 0$) is $\frac{1}{c}R_i \to R_i$.
    - The inverse of row replacement $R_i + cR_j \to R_i$ is $R_i - cR_j \to R_i$.

    Therefore, any solution of the original system satisfies all equations of the new system, and vice versa. Hence the two systems have the same solution set. $\blacksquare$

---

## 1.3 Gaussian Elimination and Gauss-Jordan Elimination

<div class="context-flow" markdown>

**From row operations to algorithms**: Elementary row operations preserve solution sets (Theorem 1.1) → systematically eliminating entries below pivots is **Gaussian elimination**; further eliminating entries above pivots is **Gauss-Jordan elimination**

</div>

### Gaussian Elimination

**Gaussian elimination** is a systematic algorithm for solving systems of linear equations. Its basic idea is to use elementary row operations to reduce the augmented matrix to row echelon form, then solve by back substitution.

**Algorithm steps**:

1. Write out the augmented matrix $[A \mid \mathbf{b}]$.
2. Starting from the leftmost column, find the first nonzero entry in that column (if the entire column is zero, move to the next column). If necessary, perform a row interchange to bring the nonzero entry to the current working row. This nonzero entry is called the **pivot**.
3. Use row replacement operations to eliminate all entries below the pivot in that column to zero.
4. Repeat steps 2–3 for the next row and next column until all rows have been processed.
5. After obtaining the row echelon form, start from the last row containing a pivot and perform **back substitution** to find each unknown.

!!! example "Example 1.3"
    Solve the system using Gaussian elimination:

    $$
    \begin{cases} x_1 + 2x_2 + x_3 = 3 \\ 2x_1 + 5x_2 + 2x_3 = 7 \\ x_1 + 3x_2 + 3x_3 = 5 \end{cases}
    $$

    Write the augmented matrix and perform row operations:

    $$
    \left(\begin{array}{ccc|c}
    1 & 2 & 1 & 3 \\
    2 & 5 & 2 & 7 \\
    1 & 3 & 3 & 5
    \end{array}\right)
    \xrightarrow{R_2 - 2R_1 \to R_2}
    \left(\begin{array}{ccc|c}
    1 & 2 & 1 & 3 \\
    0 & 1 & 0 & 1 \\
    1 & 3 & 3 & 5
    \end{array}\right)
    $$

    $$
    \xrightarrow{R_3 - R_1 \to R_3}
    \left(\begin{array}{ccc|c}
    1 & 2 & 1 & 3 \\
    0 & 1 & 0 & 1 \\
    0 & 1 & 2 & 2
    \end{array}\right)
    \xrightarrow{R_3 - R_2 \to R_3}
    \left(\begin{array}{ccc|c}
    1 & 2 & 1 & 3 \\
    0 & 1 & 0 & 1 \\
    0 & 0 & 2 & 1
    \end{array}\right)
    $$

    Back substitution: from the third row $2x_3 = 1$, we get $x_3 = \frac{1}{2}$; from the second row $x_2 = 1$; from the first row $x_1 + 2(1) + \frac{1}{2} = 3$, we get $x_1 = \frac{1}{2}$.

    The solution is $x_1 = \frac{1}{2},\; x_2 = 1,\; x_3 = \frac{1}{2}$.

### Gauss-Jordan Elimination

**Gauss-Jordan elimination** goes further than Gaussian elimination: it eliminates not only entries below each pivot but also entries above, and scales each pivot to $1$, ultimately reducing the augmented matrix to **reduced row echelon form**.

!!! example "Example 1.4"
    Continuing from the result of Example 1.3, perform Gauss-Jordan elimination:

    $$
    \left(\begin{array}{ccc|c}
    1 & 2 & 1 & 3 \\
    0 & 1 & 0 & 1 \\
    0 & 0 & 2 & 1
    \end{array}\right)
    \xrightarrow{\frac{1}{2}R_3 \to R_3}
    \left(\begin{array}{ccc|c}
    1 & 2 & 1 & 3 \\
    0 & 1 & 0 & 1 \\
    0 & 0 & 1 & \frac{1}{2}
    \end{array}\right)
    $$

    $$
    \xrightarrow{R_1 - R_3 \to R_1}
    \left(\begin{array}{ccc|c}
    1 & 2 & 0 & \frac{5}{2} \\
    0 & 1 & 0 & 1 \\
    0 & 0 & 1 & \frac{1}{2}
    \end{array}\right)
    \xrightarrow{R_1 - 2R_2 \to R_1}
    \left(\begin{array}{ccc|c}
    1 & 0 & 0 & \frac{1}{2} \\
    0 & 1 & 0 & 1 \\
    0 & 0 & 1 & \frac{1}{2}
    \end{array}\right)
    $$

    The solution $x_1 = \frac{1}{2},\; x_2 = 1,\; x_3 = \frac{1}{2}$ can be read off directly, without back substitution.

---

## 1.4 Row Echelon Form and Reduced Row Echelon Form

<div class="context-flow" markdown>

**Why standard forms are needed**: The goal of elimination needs a precise definition → **REF** suffices for back substitution, **RREF** allows reading off solutions directly → the uniqueness of RREF (Theorem 1.2) makes it the "canonical representative" of a matrix

</div>

!!! definition "Definition 1.6 (Row echelon form)"
    A matrix is in **row echelon form** (REF) if and only if it satisfies the following conditions:

    1. All zero rows are at the bottom of the matrix.
    2. The first nonzero entry of each nonzero row (called the **pivot** or **leading entry**) is strictly to the right of the pivot of the row above.

!!! definition "Definition 1.7 (Reduced row echelon form)"
    A matrix is in **reduced row echelon form** (RREF) if and only if it satisfies the following conditions:

    1. It is in row echelon form.
    2. Each pivot equals $1$ (called a **leading 1**).
    3. All other entries in each pivot column are $0$.

!!! example "Example 1.5"
    The following matrix is in row echelon form but not reduced row echelon form:

    $$
    \begin{pmatrix} 2 & 1 & 3 & 4 \\ 0 & 0 & 3 & 1 \\ 0 & 0 & 0 & 0 \end{pmatrix}
    $$

    The following matrices are in reduced row echelon form:

    $$
    \begin{pmatrix} 1 & 0 & 0 & 2 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & 3 \end{pmatrix}, \qquad
    \begin{pmatrix} 1 & 3 & 0 & 5 \\ 0 & 0 & 1 & 2 \\ 0 & 0 & 0 & 0 \end{pmatrix}
    $$

!!! theorem "Theorem 1.2 (Uniqueness of reduced row echelon form)"
    Every matrix is row equivalent to a unique reduced row echelon form matrix.

??? proof "Proof"
    **Existence**: By Gauss-Jordan elimination, any matrix can be reduced to reduced row echelon form through finitely many elementary row operations.

    **Uniqueness**: Suppose matrix $A$ is row equivalent to two reduced row echelon form matrices $R_1$ and $R_2$. Since row-equivalent matrices correspond to systems with the same solutions, $R_1$ and $R_2$ correspond to homogeneous systems with exactly the same solution set.

    We prove $R_1 = R_2$. First, $R_1$ and $R_2$ have pivot columns in the same positions (since pivot columns correspond precisely to basic variables, which are uniquely determined by the solution set). Second, for each pivot column, the reduced row echelon form has only a $1$ at the pivot position and $0$s elsewhere, so the corresponding columns are identical. For free variable columns, their values are uniquely determined by the solution set of the system. Combining these observations gives $R_1 = R_2$. $\blacksquare$

!!! note "Note"
    Row echelon form is **not** unique. The same matrix can be reduced to different row echelon forms through different sequences of elementary row operations. However, the reduced row echelon form is unique — this is an important property.

---

## 1.5 Existence and Uniqueness of Solutions

<div class="context-flow" markdown>

**From algorithm to theory**: The pivot positions in RREF completely determine the fate of solutions → **pivot column = basic variable**, non-pivot column = free variable → three outcomes: no solution / unique solution / infinitely many solutions

</div>

!!! definition "Definition 1.8 (Consistent and inconsistent)"
    If a system of linear equations has a solution (at least one), it is called **consistent**; if it has no solution, it is called **inconsistent**.

!!! definition "Definition 1.9 (Pivot position and pivot column)"
    A position in a matrix that corresponds to a pivot in a row echelon form is called a **pivot position**. A column containing a pivot position is called a **pivot column**. The unknowns corresponding to pivot columns are called **basic variables**, and the remaining unknowns are called **free variables**.

!!! theorem "Theorem 1.3 (Existence of solutions)"
    The system of linear equations $A\mathbf{x} = \mathbf{b}$ is consistent if and only if the column of $\mathbf{b}$ in the row echelon form of the augmented matrix $[A \mid \mathbf{b}]$ is not a pivot column.

    Equivalently, the system is consistent if and only if the row echelon form of the augmented matrix contains no row of the form $(0\; 0\; \cdots\; 0 \mid c)$ where $c \neq 0$.

??? proof "Proof"
    If the row echelon form contains a row $(0\; 0\; \cdots\; 0 \mid c)$ ($c \ne 0$), the corresponding equation is $0 = c$, a contradiction, so the system has no solution.

    If no such row appears, then every nonzero row has at least one pivot in a coefficient column, and by back substitution (assigning arbitrary values to free variables), at least one solution can be found, so the system is consistent. $\blacksquare$

!!! theorem "Theorem 1.4 (Uniqueness of solutions)"
    If a consistent system of linear equations has no free variables (i.e., every column is a pivot column), then the system has a unique solution; if free variables exist, then the system has infinitely many solutions.

??? proof "Proof"
    If there are no free variables, every unknown is a basic variable in the reduced row echelon form with only one $1$ in its corresponding column, so the solution is uniquely determined.

    If free variables exist, assigning different real values to the free variable yields different solutions, so there are infinitely many solutions. Since the free variables can take any value in $\mathbb{R}$, the solution set is infinite. $\blacksquare$

<div class="context-flow" markdown>

**Key insight**: A system of linear equations cannot have "exactly 2 solutions" — having 2 implies infinitely many, because the solution set is closed under affine combinations

</div>

!!! proposition "Proposition 1.1"
    The solution of a system of linear equations falls into exactly one of three cases: **no solution**, **a unique solution**, or **infinitely many solutions**. There is no system of linear equations with exactly finitely many (more than 1) solutions.

??? proof "Proof"
    Suppose the system has two distinct solutions $\mathbf{x}_1$ and $\mathbf{x}_2$. For any real number $t$, let $\mathbf{x}(t) = (1-t)\mathbf{x}_1 + t\mathbf{x}_2$. Since

    $$
    A\mathbf{x}(t) = (1-t)A\mathbf{x}_1 + tA\mathbf{x}_2 = (1-t)\mathbf{b} + t\mathbf{b} = \mathbf{b},
    $$

    $\mathbf{x}(t)$ is also a solution of the system. As $t$ ranges over all real numbers, we obtain infinitely many distinct solutions. $\blacksquare$

!!! example "Example 1.6"
    Determine the solution situation for the following system:

    $$
    \begin{cases} x_1 + 2x_2 - x_3 = 1 \\ 2x_1 + 4x_2 - 2x_3 = 3 \end{cases}
    $$

    The augmented matrix is

    $$
    \left(\begin{array}{ccc|c} 1 & 2 & -1 & 1 \\ 2 & 4 & -2 & 3 \end{array}\right) \xrightarrow{R_2 - 2R_1 \to R_2} \left(\begin{array}{ccc|c} 1 & 2 & -1 & 1 \\ 0 & 0 & 0 & 1 \end{array}\right)
    $$

    The second row is $(0\; 0\; 0 \mid 1)$, corresponding to the equation $0 = 1$, a contradiction. The system has **no solution**.

!!! example "Example 1.7"
    Solve the system

    $$
    \begin{cases} x_1 - 2x_2 + x_3 + x_4 = 3 \\ 2x_1 - 4x_2 + 3x_4 = 2 \\ -x_1 + 2x_2 - x_3 + 2x_4 = -1 \end{cases}
    $$

    Reduce the augmented matrix to reduced row echelon form:

    $$
    \left(\begin{array}{cccc|c} 1 & -2 & 1 & 1 & 3 \\ 2 & -4 & 0 & 3 & 2 \\ -1 & 2 & -1 & 2 & -1 \end{array}\right) \xrightarrow{\text{row operations}} \left(\begin{array}{cccc|c} 1 & -2 & 0 & 0 & 5 \\ 0 & 0 & 1 & 0 & -2 \\ 0 & 0 & 0 & 1 & -2 \end{array}\right)
    $$

    The pivot columns are columns 1, 3, and 4, corresponding to basic variables $x_1, x_3, x_4$. Column 2 is a non-pivot column, so $x_2$ is a free variable. Setting $x_2 = t$,

    $$
    x_1 = 5 + 2t, \quad x_2 = t, \quad x_3 = -2, \quad x_4 = -2.
    $$

---

## 1.6 Homogeneous Systems of Linear Equations

<div class="context-flow" markdown>

**From special to general**: $A\mathbf{x}=\mathbf{0}$ always has the trivial solution → focus on the existence of nontrivial solutions → the solution set is closed under addition and scalar multiplication → forms a **subspace** (→ Chapter 4 null space, Chapter 5 kernel)

</div>

!!! definition "Definition 1.10 (Homogeneous system of linear equations)"
    If all constant terms $b_i = 0$ in a system of linear equations, i.e., the system is

    $$
    A\mathbf{x} = \mathbf{0},
    $$

    it is called a **homogeneous system of linear equations**. Otherwise it is called a **nonhomogeneous system**.

A homogeneous system always has a solution $\mathbf{x} = \mathbf{0}$, called the **trivial solution**. A nonzero solution is called a **nontrivial solution**.

!!! theorem "Theorem 1.5 (Condition for nontrivial solutions)"
    The homogeneous system $A\mathbf{x} = \mathbf{0}$ has a nontrivial solution if and only if the system has free variables.

    In particular, if the number of equations $m$ is less than the number of unknowns $n$ (i.e., $m < n$), then the homogeneous system must have a nontrivial solution.

??? proof "Proof"
    A homogeneous system is always consistent ($\mathbf{x} = \mathbf{0}$ is a solution). Therefore, by Theorem 1.4, it has nontrivial solutions if and only if free variables exist.

    For the case $m < n$: the augmented matrix is an $m \times (n+1)$ matrix, and the coefficient matrix is an $m \times n$ matrix. There are at most $m$ pivots (at most one per row), so there are at least $n - m > 0$ free variables, hence nontrivial solutions exist. $\blacksquare$

!!! theorem "Theorem 1.6 (Properties of the solution space of homogeneous systems)"
    The solution set $S$ of the homogeneous system $A\mathbf{x} = \mathbf{0}$ satisfies:

    1. $\mathbf{0} \in S$.
    2. If $\mathbf{x}_1, \mathbf{x}_2 \in S$, then $\mathbf{x}_1 + \mathbf{x}_2 \in S$ (closed under addition).
    3. If $\mathbf{x}_1 \in S$, $c \in \mathbb{R}$, then $c\mathbf{x}_1 \in S$ (closed under scalar multiplication).

    Therefore $S$ is a **subspace** of $\mathbb{R}^n$, called the **null space** or **kernel** of $A$.

??? proof "Proof"
    1. $A\mathbf{0} = \mathbf{0}$, so $\mathbf{0} \in S$.

    2. If $A\mathbf{x}_1 = \mathbf{0}$ and $A\mathbf{x}_2 = \mathbf{0}$, then $A(\mathbf{x}_1 + \mathbf{x}_2) = A\mathbf{x}_1 + A\mathbf{x}_2 = \mathbf{0} + \mathbf{0} = \mathbf{0}$.

    3. If $A\mathbf{x}_1 = \mathbf{0}$, then $A(c\mathbf{x}_1) = cA\mathbf{x}_1 = c\mathbf{0} = \mathbf{0}$. $\blacksquare$

!!! example "Example 1.8"
    Find the general solution of the homogeneous system:

    $$
    \begin{cases} x_1 + 2x_2 - x_3 + x_4 = 0 \\ 2x_1 + 4x_2 + x_3 - 2x_4 = 0 \\ 3x_1 + 6x_2 + 2x_4 = 0 \end{cases}
    $$

    Reduce the augmented matrix to reduced row echelon form:

    $$
    \left(\begin{array}{cccc|c} 1 & 2 & -1 & 1 & 0 \\ 2 & 4 & 1 & -2 & 0 \\ 3 & 6 & 0 & 2 & 0 \end{array}\right)
    \xrightarrow{\text{row operations}}
    \left(\begin{array}{cccc|c} 1 & 2 & 0 & -\frac{1}{3} & 0 \\ 0 & 0 & 1 & -\frac{4}{3} & 0 \\ 0 & 0 & 0 & 0 & 0 \end{array}\right)
    $$

    The pivot columns are columns 1 and 3. Free variables: $x_2 = s$, $x_4 = t$. General solution:

    $$
    \begin{pmatrix} x_1 \\ x_2 \\ x_3 \\ x_4 \end{pmatrix}
    = s \begin{pmatrix} -2 \\ 1 \\ 0 \\ 0 \end{pmatrix}
    + t \begin{pmatrix} \frac{1}{3} \\ 0 \\ \frac{4}{3} \\ 1 \end{pmatrix}, \quad s, t \in \mathbb{R}.
    $$

---

## 1.7 Structure of Solutions

<div class="context-flow" markdown>

**Unified perspective**: General solution of nonhomogeneous system = particular solution $\mathbf{x}_0$ + general solution of homogeneous system → the solution set is a **translate** of a subspace (affine subspace) → this is the prototype of the superposition principle

</div>

!!! theorem "Theorem 1.7 (Structure of solutions of nonhomogeneous systems)"
    Let $A\mathbf{x} = \mathbf{b}$ be a consistent nonhomogeneous system and $\mathbf{x}_0$ be one of its **particular solutions**. Then the general solution of the system can be expressed as

    $$
    \mathbf{x} = \mathbf{x}_0 + \mathbf{x}_h,
    $$

    where $\mathbf{x}_h$ is the general solution of the corresponding homogeneous system $A\mathbf{x} = \mathbf{0}$.

    In other words, **the general solution of a nonhomogeneous system = a particular solution + the general solution of the homogeneous system**.

??? proof "Proof"
    **(Every solution has this form)** Let $\mathbf{x}_1$ be any solution of $A\mathbf{x} = \mathbf{b}$. Then

    $$
    A(\mathbf{x}_1 - \mathbf{x}_0) = A\mathbf{x}_1 - A\mathbf{x}_0 = \mathbf{b} - \mathbf{b} = \mathbf{0},
    $$

    so $\mathbf{x}_h = \mathbf{x}_1 - \mathbf{x}_0$ is a solution of the homogeneous system, and thus $\mathbf{x}_1 = \mathbf{x}_0 + \mathbf{x}_h$.

    **(Every expression of this form is a solution)** Conversely, if $A\mathbf{x}_h = \mathbf{0}$, then

    $$
    A(\mathbf{x}_0 + \mathbf{x}_h) = A\mathbf{x}_0 + A\mathbf{x}_h = \mathbf{b} + \mathbf{0} = \mathbf{b},
    $$

    so $\mathbf{x}_0 + \mathbf{x}_h$ is indeed a solution of the nonhomogeneous system. $\blacksquare$

!!! note "Note"
    This theorem reveals the geometric structure of the solution set of a linear system. The solution set of the homogeneous system $A\mathbf{x} = \mathbf{0}$ is a subspace passing through the origin (a line, plane, etc.), while the solution set of the nonhomogeneous system $A\mathbf{x} = \mathbf{b}$ is a **translate** of this subspace (an affine subspace), i.e., the subspace shifted to the position of the particular solution $\mathbf{x}_0$.

!!! example "Example 1.9"
    Solve the system

    $$
    \begin{cases} x_1 + 2x_2 + x_3 = 4 \\ 2x_1 + 4x_2 + 3x_3 = 9 \end{cases}
    $$

    Reduce the augmented matrix:

    $$
    \left(\begin{array}{ccc|c} 1 & 2 & 1 & 4 \\ 2 & 4 & 3 & 9 \end{array}\right)
    \xrightarrow{R_2 - 2R_1}
    \left(\begin{array}{ccc|c} 1 & 2 & 1 & 4 \\ 0 & 0 & 1 & 1 \end{array}\right)
    \xrightarrow{R_1 - R_2}
    \left(\begin{array}{ccc|c} 1 & 2 & 0 & 3 \\ 0 & 0 & 1 & 1 \end{array}\right)
    $$

    **Particular solution** (set free variable $x_2 = 0$): $\mathbf{x}_0 = (3, 0, 1)^T$.

    **Homogeneous general solution** ($x_2 = t$): $\mathbf{x}_h = t(-2, 1, 0)^T$.

    **General solution**:

    $$
    \mathbf{x} = \begin{pmatrix} 3 \\ 0 \\ 1 \end{pmatrix} + t \begin{pmatrix} -2 \\ 1 \\ 0 \end{pmatrix}, \quad t \in \mathbb{R}.
    $$

!!! example "Example 1.10"
    Suppose $A$ is a $4 \times 5$ matrix with rank $3$, and $A\mathbf{x} = \mathbf{b}$ has a solution. Analyze the structure of the solutions.

    Since $A$ has $5$ unknowns and rank $3$ (i.e., $3$ pivot columns), there are $5 - 3 = 2$ free variables. The solution space of $A\mathbf{x} = \mathbf{0}$ is a $2$-dimensional subspace of $\mathbb{R}^5$; let its fundamental system of solutions be $\{\mathbf{x}_1, \mathbf{x}_2\}$.

    The general solution of the nonhomogeneous system is

    $$
    \mathbf{x} = \mathbf{x}_0 + c_1 \mathbf{x}_1 + c_2 \mathbf{x}_2, \quad c_1, c_2 \in \mathbb{R},
    $$

    where $\mathbf{x}_0$ is a particular solution. Geometrically, this is a $2$-dimensional affine subspace in $\mathbb{R}^5$ (a generalization of a plane).
