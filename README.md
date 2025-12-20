### Project Overview
In this project, various important methods of Numerical Methods have been implemented and presented in a topic-wise manner. The project is divided into major sections including Solution of Linear Equations, Non-Linear Equations, Differential Equation Solving, Interpolation Methods, Numerical Differentiation, Curve Fitting/Regression, and Numerical Integration. Each method under every section is organized separately into Theory, C++ Code, Input, and Output parts. The overall structure and topic-wise organization of the project are clearly documented in the README.md file, making the work easy to understand and follow. This project serves as a well-organized and user-friendly reference for learning numerical methods. 
---

# Table of Contents

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)
      
- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Newton Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)

- [Differential Equation Solving](#differential-equation-solving)
  - [Runge-Kutta 4th Order Method](#runge-kutta-4th-order-method)
    - [Theory](#runge-kutta-theory)
    - [Code](#runge-kutta-code)
    - [Input](#runge-kutta-input)
    - [Output](#runge-kutta-output)

- [Interpolation Methods](#interpolation-methods)
  - [Newton Forward Interpolation](#newton-forward-interpolation)
    - [Theory](#newton-forward-theory)
    - [Code](#newton-forward-code)
    - [Input](#newton-forward-input)
    - [Output](#newton-forward-output)
  - [Newton Backward Interpolation](#newton-backward-interpolation)
    - [Theory](#newton-backward-theory)
    - [Code](#newton-backward-code)
    - [Input](#newton-backward-input)
    - [Output](#newton-backward-output)
  - [Newton Divided Difference Interpolation](#newton-divided-difference-interpolation)
    - [Theory](#newton-divided-difference-theory)
    - [Code](#newton-divided-difference-code)
    - [Input](#newton-divided-difference-input)
    - [Output](#newton-divided-difference-output)

- [Numerical Differentiation](#numerical-differentiation)
  - [Differentiation by Forward Interpolation](#differentiation-by-forward-interpolation)
    - [Theory](#differentiation-forward-theory)
    - [Code](#differentiation-forward-code)
    - [Input](#differentiation-forward-input)
    - [Output](#differentiation-forward-output)
  - [Differentiation by Backward Interpolation](#differentiation-by-backward-interpolation)
    - [Theory](#differentiation-backward-theory)
    - [Code](#differentiation-backward-code)
    - [Input](#differentiation-backward-input)
    - [Output](#differentiation-backward-output)

- [Curve Fitting / Regression](#curve-fitting--regression)
  - [Linear Regression](#linear-regression)
    - [Theory](#linear-regression-theory)
    - [Code](#linear-regression-code)
    - [Input](#linear-regression-input)
    - [Output](#linear-regression-output)
  - [Polynomial Regression](#polynomial-regression)
    - [Theory](#polynomial-regression-theory)
    - [Code](#polynomial-regression-code)
    - [Input](#polynomial-regression-input)
    - [Output](#polynomial-regression-output)
  - [Transcendental Regression](#transcendental-regression)
    - [Theory](#transcendental-regression-theory)
    - [Code](#transcendental-regression-code)
    - [Input](#transcendental-regression-input)
    - [Output](#transcendental-regression-output)

- [Numerical Integration](#numerical-integration)
  - [Simpson's 1/3 Rule](#simpsons-13-rule)
    - [Theory](#simpson-13-theory)
    - [Code](#simpson-13-code)
    - [Input](#simpson-13-input)
    - [Output](#simpson-13-output)
  - [Simpson's 3/8 Rule](#simpsons-38-rule)
    - [Theory](#simpson-38-theory)
    - [Code](#simpson-38-code)
    - [Input](#simpson-38-input)
    - [Output](#simpson-38-output)

---

# Numerical Methods

---


## Solution of Linear Equations

---

### Gauss Elimination Method

#### Gauss Elimination Theory

### Algorithm
1. Write the system of equations in augmented matrix form.
2. Use row operations to eliminate variables below the main diagonal.
3. Convert the matrix into upper triangular form.
4. Continue elimination until the last equation contains one unknown.
5. Apply backward substitution to find the solution.

### Mathematical Formula
For a system AX = B, the augmented matrix is reduced as:

[A | B] → [U | C]

where U is an upper triangular matrix.

Backward substitution:

x_n = c_n / u_nn

x_i = ( c_i − ( u_i(i+1)x_(i+1) + u_i(i+2)x_(i+2) + … + u_in x_n ) ) / u_ii

### Brief
Gauss Elimination is a basic numerical technique for solving linear systems.
It works by simplifying equations step by step using row operations.
The method follows a clear and logical structure.
It is widely used in engineering and science problems.
Care must be taken when pivot elements are very small.


#### Gauss Elimination Code

```cpp
#include<bits/stdc++.h>
using namespace std;

void gaussElimination(ifstream &in,ofstream &out)
{
    for(int c=1;c<=3;c++)
    {
        int n;
        in>>n;

        vector<vector<double>> a(n,vector<double>(n+1));
        for(int i=0;i<n;i++)
            for(int j=0;j<=n;j++)
                in>>a[i][j];

        for(int i=0;i<n;i++)
        {
            int p=i;
            for(int j=i;j<n;j++)
                if(fabs(a[j][i])>fabs(a[p][i])) p=j;
            swap(a[i],a[p]);

            if(fabs(a[i][i])<1e-9) continue;

            for(int j=i+1;j<n;j++)
            {
                double r=a[j][i]/a[i][i];
                for(int k=i;k<=n;k++)
                    a[j][k]-=r*a[i][k];
            }
        }

        int r1=0,r2=0;
        for(int i=0;i<n;i++)
        {
            bool nz=false;
            for(int j=0;j<n;j++)
                if(fabs(a[i][j])>1e-9) nz=true;
            if(nz) r1++;

            nz=false;
            for(int j=0;j<=n;j++)
                if(fabs(a[i][j])>1e-9) nz=true;
            if(nz) r2++;
        }

        out<<"Case "<<c<<": ";

        if(r1!=r2)
        {
            out<<"No Solution\n";
        }
        else if(r1<n)
        {
            out<<"Infinite Solutions\n";
        }
        else
        {
            vector<double> x(n);
            for(int i=n-1;i>=0;i--)
            {
                x[i]=a[i][n];
                for(int j=i+1;j<n;j++)
                    x[i]-=a[i][j]*x[j];
                x[i]/=a[i][i];
            }

            out<<"Unique Solution ";
            for(int i=0;i<n;i++)
                out<<fixed<<setprecision(3)<<x[i]<<" ";
            out<<"\n";
        }
    }
}

int main()
{
    ifstream in("gauss_input.txt");
    ofstream out("gauss_output.txt");
    gaussElimination(in,out);
}

```

#### Gauss Elimination Input

```
3
2
1 1 2
2 2 5
2
1 1 2
2 2 4
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3

```

#### Gauss Elimination Output

```
Case 1: No Solution
Case 2: Infinite Solutions
Case 3: Unique Solution 2.000 3.000 -1.000
```

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory

### Algorithm
1. Convert the system into augmented matrix form.
2. Make the leading element of each row equal to one.
3. Eliminate all other elements in the pivot column.
4. Reduce the matrix completely.
5. Read the solution directly.

### Mathematical Formula

[A | B] → [I | X]

where I is the identity matrix and X is the solution vector.

### Brief
Gauss Jordan Elimination is an extension of Gauss Elimination.
It removes the need for backward substitution.
Each variable is isolated directly during elimination.
The solution is easy to interpret.
Although slower, it is very systematic and accurate.


#### Gauss Jordan Code

```cpp
#include<bits/stdc++.h>
using namespace std;

void gaussJordan(ifstream &in,ofstream &out)
{
    for(int c=1;c<=3;c++)
    {
        int n;
        in>>n;

        vector<vector<double>> a(n,vector<double>(n+1));
        for(int i=0;i<n;i++)
            for(int j=0;j<=n;j++)
                in>>a[i][j];

        for(int i=0;i<n;i++)
        {
            int p=i;
            for(int j=i;j<n;j++)
                if(fabs(a[j][i])>fabs(a[p][i])) p=j;
            swap(a[i],a[p]);

            if(fabs(a[i][i])<1e-9) continue;

            double d=a[i][i];
            for(int j=0;j<=n;j++)
                a[i][j]/=d;

            for(int k=0;k<n;k++)
                if(k!=i)
                {
                    double r=a[k][i];
                    for(int j=0;j<=n;j++)
                        a[k][j]-=r*a[i][j];
                }
        }

        out<<"Case "<<c<<":\n";
        out<<"Reduced Row Echelon Form:\n";
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<=n;j++)
                out<<fixed<<setprecision(3)<<a[i][j]<<" ";
            out<<"\n";
        }

        int r1=0,r2=0;
        for(int i=0;i<n;i++)
        {
            bool nz=false;
            for(int j=0;j<n;j++)
                if(fabs(a[i][j])>1e-9) nz=true;
            if(nz) r1++;

            nz=false;
            for(int j=0;j<=n;j++)
                if(fabs(a[i][j])>1e-9) nz=true;
            if(nz) r2++;
        }

        if(r1!=r2)
        {
            out<<"No Solution\n\n";
        }
        else if(r1<n)
        {
            out<<"Infinite Solutions\n\n";
        }
        else
        {
            out<<"Unique Solution ";
            for(int i=0;i<n;i++)
                out<<fixed<<setprecision(3)<<a[i][n]<<" ";
            out<<"\n\n";
        }
    }
}

int main()
{
    ifstream in("gj_input.txt");
    ofstream out("gj_output.txt");
    gaussJordan(in,out);
}

```

#### Gauss Jordan Input

```
3
2
1 1 2
2 2 5
2
1 1 2
2 2 4
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3

```

#### Gauss Jordan Output

```
Case 1:
Reduced Row Echelon Form:
1.000 1.000 0.000
0.000 0.000 1.000
No Solution

Case 2:
Reduced Row Echelon Form:
1.000 1.000 2.000
0.000 0.000 0.000
Infinite Solutions

Case 3:
Reduced Row Echelon Form:
1.000 0.000 0.000 2.000
0.000 1.000 0.000 3.000
0.000 0.000 1.000 -1.000
Unique Solution 2.000 3.000 -1.000

```

---

### LU Decomposition Method

#### LU Decomposition Theory

### Algorithm
1. Decompose the matrix A into L and U.
2. Solve LY = B using forward substitution.
3. Solve UX = Y using backward substitution.
4. Obtain the final solution.

### Mathematical Formula

A = L × U

Forward substitution:

y_1 = b_1
y_i = b_i − ( l_i1 y_1 + l_i2 y_2 + … + l_i(i−1) y_(i−1) )

Backward substitution:

x_n = y_n / u_nn

x_i = ( y_i − ( u_i(i+1)x_(i+1) + … + u_in x_n ) ) / u_ii

### Brief
LU Decomposition simplifies solving large systems.
It reduces repeated computations.
The method is efficient for multiple right-hand sides.
It is commonly used in numerical software.
LU decomposition improves computational speed.


#### LU Decomposition Code

```cpp
#include<bits/stdc++.h>
using namespace std;

void luDecomposition(ifstream &in,ofstream &out)
{
    for(int c=1;c<=3;c++)
    {
        int n;
        in>>n;

        vector<vector<double>> a(n,vector<double>(n));
        vector<double> b(n);

        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
                in>>a[i][j];
        for(int i=0;i<n;i++)
            in>>b[i];

        vector<vector<double>> l(n,vector<double>(n,0)),u(n,vector<double>(n,0));

        for(int i=0;i<n;i++)
        {
            for(int k=i;k<n;k++)
            {
                double s=0;
                for(int j=0;j<i;j++)
                    s+=l[i][j]*u[j][k];
                u[i][k]=a[i][k]-s;
            }

            for(int k=i;k<n;k++)
            {
                if(i==k) l[i][i]=1;
                else
                {
                    if(fabs(u[i][i])<1e-9) continue;
                    double s=0;
                    for(int j=0;j<i;j++)
                        s+=l[k][j]*u[j][i];
                    l[k][i]=(a[k][i]-s)/u[i][i];
                }
            }
        }

        int r1=0,r2=0;
        for(int i=0;i<n;i++)
        {
            bool nz=false;
            for(int j=0;j<n;j++)
                if(fabs(a[i][j])>1e-9) nz=true;
            if(nz) r1++;

            if(nz || fabs(b[i])>1e-9) r2++;
        }

        out<<"Case "<<c<<":\n";

        out<<"Lower Triangular Matrix (L):\n";
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
                out<<fixed<<setprecision(3)<<l[i][j]<<" ";
            out<<"\n";
        }

        out<<"Upper Triangular Matrix (U):\n";
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
                out<<fixed<<setprecision(3)<<u[i][j]<<" ";
            out<<"\n";
        }

        if(r1!=r2)
        {
            out<<"No Solution\n\n";
            continue;
        }
        else if(r1<n)
        {
            out<<"Infinite Solutions\n\n";
            continue;
        }

        vector<double> y(n),x(n);

        for(int i=0;i<n;i++)
        {
            y[i]=b[i];
            for(int j=0;j<i;j++)
                y[i]-=l[i][j]*y[j];
        }

        for(int i=n-1;i>=0;i--)
        {
            x[i]=y[i];
            for(int j=i+1;j<n;j++)
                x[i]-=u[i][j]*x[j];
            x[i]/=u[i][i];
        }

        out<<"Unique Solution ";
        for(int i=0;i<n;i++)
            out<<fixed<<setprecision(3)<<x[i]<<" ";
        out<<"\n\n";
    }
}

int main()
{
    ifstream in("lu_input.txt");
    ofstream out("lu_output.txt");
    luDecomposition(in,out);
}

```

#### LU Decomposition Input

```
3
2
1 1
2 2
2 5
3
1 1 1
2 2 2
3 3 3
6 12 20
3
2 1 1
4 -6 0
-2 7 2
5 -2 9

```

#### LU Decomposition Output

```
Case 1:
Lower Triangular Matrix (L):
1.000 0.000
2.000 1.000
Upper Triangular Matrix (U):
1.000 1.000
0.000 0.000
No Solution

Case 2:
Lower Triangular Matrix (L):
1.000 0.000 0.000
2.000 1.000 0.000
3.000 0.000 1.000
Upper Triangular Matrix (U):
1.000 1.000 1.000
0.000 0.000 0.000
0.000 0.000 0.000
Infinite Solutions

Case 3:
Lower Triangular Matrix (L):
1.000 0.000 0.000
2.000 1.000 0.000
-1.000 -1.000 1.000
Upper Triangular Matrix (U):
2.000 1.000 1.000
0.000 -8.000 -2.000
0.000 0.000 1.000
Unique Solution 1.000 1.000 2.000

```

---

### Matrix Inversion

#### Matrix Inversion Theory

### Algorithm
1. Check whether the determinant of A is non-zero.
2. Form the augmented matrix [A | I].
3. Apply Gauss Jordan elimination.
4. Extract the inverse matrix.
5. Multiply inverse with B to get the solution.

### Mathematical Formula

X = A_inverse × B

Matrix transformation:

[A | I] → [I | A_inverse]

### Brief
Matrix Inversion Method gives a direct solution.
It is easy to understand mathematically.
The method works only for non-singular matrices.
It becomes inefficient for large matrices.
Inverse computation is computationally expensive.


#### Matrix Inversion Code

```cpp
#include<bits/stdc++.h>
using namespace std;

void matrixInversion(ifstream &in,ofstream &out)
{
    for(int c=1;c<=3;c++)
    {
        int n;
        in>>n;

        vector<vector<double>> a(n,vector<double>(2*n));
        vector<double> b(n);

        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
                in>>a[i][j];

        for(int i=0;i<n;i++)
            in>>b[i];

        for(int i=0;i<n;i++)
            for(int j=n;j<2*n;j++)
                a[i][j]=(i==j-n);

        out<<"Case "<<c<<":\n";

        out<<"Vector B:\n";
        for(int i=0;i<n;i++)
            out<<fixed<<setprecision(3)<<b[i]<<" ";
        out<<"\n";

        out<<"Augmented Matrix [A | I]:\n";
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<2*n;j++)
                out<<fixed<<setprecision(3)<<a[i][j]<<" ";
            out<<"\n";
        }

        for(int i=0;i<n;i++)
        {
            int p=i;
            for(int j=i;j<n;j++)
                if(fabs(a[j][i])>fabs(a[p][i])) p=j;
            swap(a[i],a[p]);

            double d=a[i][i];
            for(int j=0;j<2*n;j++)
                a[i][j]/=d;

            for(int k=0;k<n;k++)
                if(k!=i)
                {
                    double r=a[k][i];
                    for(int j=0;j<2*n;j++)
                        a[k][j]-=r*a[i][j];
                }
        }

        out<<"Reduced Row Echelon Form:\n";
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<2*n;j++)
                out<<fixed<<setprecision(3)<<a[i][j]<<" ";
            out<<"\n";
        }

        out<<"Inverse Matrix:\n";
        for(int i=0;i<n;i++)
        {
            for(int j=n;j<2*n;j++)
                out<<fixed<<setprecision(3)<<a[i][j]<<" ";
            out<<"\n";
        }

        vector<double> x(n,0);
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
                x[i]+=a[i][j+n]*b[j];

        out<<"Solution Vector:\n";
        for(int i=0;i<n;i++)
            out<<fixed<<setprecision(3)<<x[i]<<" ";
        out<<"\n\n";
    }
}

int main()
{
    ifstream in("inv_input.txt");
    ofstream out("inv_output.txt");
    matrixInversion(in,out);
}

```

#### Matrix Inversion Input

```
3
2
1 2
2 4
5 10
3
1 1 1
2 2 2
3 3 3
6 12 18
3
2 1 1
4 -6 0
-2 7 2
5 -2 9

```

#### Matrix Inversion Output

```
Case 1:
Vector B:
5.000 10.000
Augmented Matrix [A | I]:
1.000 2.000 1.000 0.000
2.000 4.000 0.000 1.000
Reduced Row Echelon Form:
1.000 2.000 0.500 -0.500
0.000 0.000 -1.000 1.000
Matrix is Singular

Case 2:
Vector B:
6.000 12.000 18.000
Augmented Matrix [A | I]:
1.000 1.000 1.000 1.000 0.000 0.000
2.000 2.000 2.000 0.000 1.000 0.000
3.000 3.000 3.000 0.000 0.000 1.000
Reduced Row Echelon Form:
1.000 1.000 1.000 0.333 -0.167 0.000
0.000 0.000 0.000 -0.667 0.333 0.000
0.000 0.000 0.000 -1.000 0.000 1.000
Matrix is Singular

Case 3:
Vector B:
5.000 -2.000 9.000
Augmented Matrix [A | I]:
2.000 1.000 1.000 1.000 0.000 0.000
4.000 -6.000 0.000 0.000 1.000 0.000
-2.000 7.000 2.000 0.000 0.000 1.000
Reduced Row Echelon Form:
1.000 0.000 0.000 0.250 0.063 -0.063
0.000 1.000 0.000 0.500 0.125 0.125
0.000 0.000 1.000 -0.750 -0.125 0.375
Inverse Matrix:
0.250 0.063 -0.063
0.500 0.125 0.125
-0.750 -0.125 0.375
Solution Vector:
1.000 1.000 2.000

```

---

## Solution of Non-Linear Equations

---

### Bisection Method

#### Bisection Theory

The Bisection Method is a numerical method used to find the roots of a continuous function. It works by repeatedly dividing an interval `[a, b]` in half and selecting the subinterval in which the function changes sign, thereby narrowing down to the root.

 
- Requires the function `f(x)` to be **continuous** on `[a, b]`.  
- The root exists only if `f(a)` and `f(b)` have **opposite signs** (`f(a)*f(b) < 0`).  
- Iteratively reduces the interval to approximate the root.  
- Simple and reliable, but may be slower compared to other methods.

**Algorithm**  
1. Choose initial interval `[a, b]` such that `f(a)*f(b) < 0`.  
2. Compute midpoint: `c = (a + b)/2`.  
3. If `f(c) = 0` (or |f(c)| < tolerance), `c` is the root.  
4. If `f(a)*f(c) < 0`, set `b = c`; else set `a = c`.  
5. Repeat steps 2–4 until desired accuracy is achieved.

**Example:**  
Find a root of `f(x) = x^2 - 4` in `[1, 3]`:  
1. Initial interval: `[a, b] = [1, 3]`  
2. Midpoint: `c = (1 + 3)/2 = 2`  
3. Check `f(c) = 2^2 - 4 = 0`  
4. Root found: `x = 2`



#### Bisection Code
```cpp
#include<bits/stdc++.h>
using namespace std;

typedef double db;


db evalPoly(db x, vector<db> &coeff){
    db res = 0, p = 1;
    for(int i = 0; i < coeff.size(); i++){
        res += coeff[i]*p;
        p *= x;
    }
    return res;
}


pair<db,int> bisectionRoot(vector<db>& coeff, db a, db b, db tol, int maxIter){
    db fa = evalPoly(a, coeff);
    db fb = evalPoly(b, coeff);
    if(fa*fb > 0) return {nan(""),0}; 
    db c;
    int iter;
    for(iter=1; iter<=maxIter; iter++){
        c = (a+b)/2.0;
        db fc = evalPoly(c, coeff);
        if(fabs(fc) < tol) return {c, iter};
        if(fa*fc < 0){ b=c; fb=fc; }
        else{ a=c; fa=fc; }
    }
    return {c, iter};
}

int main(){
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int t; fin>>t;
    fout<<"Total Test Cases: "<<t<<"\n\n";

    for(int tc=1; tc<=t; tc++){
        int degree; fin>>degree;
        vector<db> coeff(degree+1);
        for(int i=0; i<=degree; i++) fin>>coeff[i];

        db a,b,tol; int maxIter;
        fin>>a>>b>>tol>>maxIter;

        fout<<"Test Case #"<<tc<<"\n";
        fout<<"Degree: "<<degree<<"\nFunction: ";
        for(int i=degree; i>=0; i--){
            if(coeff[i] == (int)coeff[i]) fout<<(int)coeff[i];
            else fout<<fixed<<setprecision(6)<<coeff[i];
            fout<<"x^"<<i;
            if(i>0) fout<<" + ";
        }
        fout<<"\nError Tolerance: "<<tol;
        fout<<"\nSearch Limit: ["<<a<<", "<<b<<"]\nRoots:\n";

        vector<pair<db,int>> roots;
        int steps=50; 
        db stepSize = (b-a)/steps;
        for(db x=a; x<b; x+=stepSize){
            auto r = bisectionRoot(coeff, x, x+stepSize, tol, maxIter);
            if(!isnan(r.first)){
                bool unique=true;
                for(auto &p: roots){
                    if(fabs(p.first - r.first) < 1e-4){ unique=false; break; }
                }
                if(unique) roots.push_back(r);
            }
        }

        for(int i=0; i<roots.size(); i++){
            fout<<"  x"<<i+1<<" = "<<fixed<<setprecision(6)<<roots[i].first<<"\n";
        }
        fout<<"\n";
    }
    fin.close(); fout.close();
    return 0;
}

````

#### Bisection Input

```
3
2
2 -3 1
0 3 0.0001 100
3
5 -6 0 1
-2 4 0.0001 100
2
3 -5 2
0 3 0.0001 100


```

#### Bisection Output

```
Total Test Cases: 3

Test Case #1
Degree: 2
Function: 1x^2 + -3x^1 + 2x^0
Error Tolerance: 0.0001
Search Limit: [0, 3]
Roots:
  x1 = 1.000078
  x2 = 1.999922

Test Case #2
Degree: 3
Function: 1x^3 + 0x^2 + -6x^1 + 5x^0
Error Tolerance: 0.000100
Search Limit: [-2.000000, 4.000000]
Roots:
  x1 = 0.999971
  x2 = 1.791309

Test Case #3
Degree: 2
Function: 2x^2 + -5x^1 + 3x^0
Error Tolerance: 0.000100
Search Limit: [0.000000, 3.000000]
Roots:
  x1 = 1.000078
  x2 = 1.499941



```

---

### False Position Method

#### False Position Theory
  
The False Position Method, also called Regula Falsi Method, is a numerical technique used to find the roots of a continuous function. It improves upon the Bisection Method by using a straight line (secant) through the interval endpoints to approximate the root, leading to faster convergence in some cases.
 
- Requires the function `f(x)` to be **continuous** on `[a, b]`.  
- The root exists only if `f(a)` and `f(b)` have **opposite signs** (`f(a)*f(b) < 0`).  
- Uses a **linear approximation** to find the next estimate of the root.  
- Convergence may be faster than the Bisection Method if the function is well-behaved.

**Algorithm:**  
1. Choose initial interval `[a, b]` such that `f(a)*f(b) < 0`.  
2. Compute the next approximation:  
c = b - (f(b)*(b - a)) / (f(b) - f(a))
3. If `f(c) = 0` (or |f(c)| < tolerance), `c` is the root.  
4. If `f(a)*f(c) < 0`, set `b = c`; else set `a = c`.  
5. Repeat steps 2–4 until desired accuracy is achieved.

**Example:**  
Find a root of `f(x) = x^2 - 4` in `[1, 3]`:  
1. Initial interval: `[a, b] = [1, 3]`  
2. Compute `c`:  
c = 3 - (f(3)(3 - 1)) / (f(3) - f(1)) = 3 - (52)/(5 - (-3)) = 3 - 10/8 = 1.75
3. Update interval based on sign of `f(c)` and repeat until convergence.



#### False Position Code

```cpp
#include<bits/stdc++.h>
using namespace std;

typedef double db;


db evalPoly(db x, vector<db> &coeff){
    db res=0, p=1;
    for(int i=0;i<coeff.size();i++){
        res+=coeff[i]*p;
        p*=x;
    }
    return res;
}


pair<db,int> falsePosition(vector<db>& coeff, db a, db b, db tol, int maxIter){
    db fa=evalPoly(a, coeff);
    db fb=evalPoly(b, coeff);
    if(fa*fb>0) return {nan(""),0};
    db c;
    int iter;
    for(iter=1; iter<=maxIter; iter++){
        c = b - fb*(b-a)/(fb-fa);
        db fc = evalPoly(c, coeff);
        if(fabs(fc)<tol) return {c, iter};
        if(fa*fc < 0){ b=c; fb=fc; }
        else{ a=c; fa=fc; }
    }
    return {c, iter};
}

int main(){
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int t; fin>>t;
    fout<<"Number of Test Cases: "<<t<<"\n\n";

    for(int tc=1; tc<=t; tc++){
        int degree; fin>>degree;
        vector<db> coeff(degree+1);
        for(int i=0;i<=degree;i++) fin>>coeff[i];

        db a,b,tol; int maxIter;
        fin>>a>>b>>tol>>maxIter;

        fout<<"--- TEST CASE "<<tc<<" ---\n";
        fout<<"Interval of search: ["<<a<<", "<<b<<"]\n";
        fout<<"Polynomial Degree: "<<degree<<"\n";
        fout<<"Given Polynomial: ";
        for(int i=degree;i>=0;i--){
            if(coeff[i]==(int)coeff[i]) fout<<(int)coeff[i];
            else fout<<fixed<<setprecision(6)<<coeff[i];
            fout<<"x^"<<i;
            if(i>0) fout<<" + ";
        }
        fout<<"\nAllowed Error: "<<tol<<"\nRoots Discovered:\n";

        vector<pair<db,int>> roots;
        int steps=50;
        db stepSize=(b-a)/steps;
        for(db x=a; x<b; x+=stepSize){
            auto r=falsePosition(coeff,x,x+stepSize,tol,maxIter);
            if(!isnan(r.first)){
                bool unique=true;
                for(auto &p:roots){
                    if(fabs(p.first-r.first)<1e-4){ unique=false; break; }
                }
                if(unique) roots.push_back(r);
            }
        }

        for(int i=0;i<roots.size();i++){
            db l = max(a, roots[i].first - stepSize/2);
            db u = min(b, roots[i].first + stepSize/2);
            fout<<"  -> Approx. Root "<<i+1<<" lies in ["<<fixed<<setprecision(3)<<l<<", "<<u<<"] = "<<fixed<<setprecision(6)<<roots[i].first<<"\n";
        }

        fout<<"Summary of All Roots:\n";
        for(int i=0;i<roots.size();i++){
            fout<<"   Root#"<<i+1<<": "<<fixed<<setprecision(6)<<roots[i].first<<"\n";
        }
        fout<<"\n";
    }

    fin.close(); fout.close();
    return 0;
}


```

#### False Position Input

```
3
2
3 -4 1
0 3 0.0001 100
3
1 -5 -3 2
-3 3 0.0001 100
2
-3 2 1
-2 2 0.0001 100


```

#### False Position Output

```
Number of Test Cases: 3

--- TEST CASE 1 ---
Interval of search: [0, 3]
Polynomial Degree: 2
Given Polynomial: 1x^2 + -4x^1 + 3x^0
Allowed Error: 0.0001
Roots Discovered:
  -> Approx. Root 1 lies in [0.970, 1.030] = 1.000008
  -> Approx. Root 2 lies in [2.970, 3.000] = 3.000000
Summary of All Roots:
   Root#1: 1.000008
   Root#2: 3.000000

--- TEST CASE 2 ---
Interval of search: [-3.000000, 3.000000]
Polynomial Degree: 3
Given Polynomial: 2x^3 + -3x^2 + -5x^1 + 1x^0
Allowed Error: 0.000100
Roots Discovered:
  -> Approx. Root 1 lies in [-1.183, -1.063] = -1.122904
  -> Approx. Root 2 lies in [0.122, 0.242] = 0.182455
  -> Approx. Root 3 lies in [2.380, 2.500] = 2.440449
Summary of All Roots:
   Root#1: -1.122904
   Root#2: 0.182455
   Root#3: 2.440449

--- TEST CASE 3 ---
Interval of search: [-2.000000, 2.000000]
Polynomial Degree: 2
Given Polynomial: 1x^2 + 2x^1 + -3x^0
Allowed Error: 0.000100
Roots Discovered:
  -> Approx. Root 1 lies in [0.960, 1.040] = 0.999996
Summary of All Roots:
   Root#1: 0.999996




```

---

### Newton Raphson Method

#### Newton Raphson Theory

The Newton-Raphson Method is a numerical technique used to find the roots of a real-valued function. It uses the tangent line at an initial guess to approximate the root and iteratively refines the estimate. It is faster than Bisection and False Position methods for well-behaved functions.
  
- Requires the function `f(x)` to be **differentiable**.  
- Convergence is **quadratic** if the initial guess is close to the root.  
- Iteratively uses the derivative of the function to refine the root estimate.  
- Sensitive to the choice of initial guess; poor guesses may lead to divergence.

**Algorithm:**  
1. Choose an initial guess `x0`.  
2. Compute the next approximation using:  
x_{n+1} = x_n - f(x_n)/f'(x_n)
3. Repeat step 2 until |f(x_n)| < tolerance or the desired accuracy is achieved.

**Example:**  
Find a root of `f(x) = x^2 - 4` starting with `x0 = 3`:  
1. Compute derivative: `f'(x) = 2x`  
2. Apply formula:  
x1 = x0 - f(x0)/f'(x0) = 3 - (3^2 - 4)/(2*3) = 3 - 5/6 = 2.1667
3. Repeat:  
x2 = 2.1667 - (2.1667^2 - 4)/(2*2.1667) ≈ 2.0064
x3 ≈ 2.0000

Root converges to `x ≈ 2`.




#### Newton Raphson Code

```cpp
#include<bits/stdc++.h>
using namespace std;
typedef double db;


db evalPoly(db x, vector<db> &coeff){
    db res=0,p=1;
    for(int i=0;i<coeff.size();i++){
        res+=coeff[i]*p;
        p*=x;
    }
    return res;
}


db evalPolyDerivative(db x, vector<db> &coeff){
    db res=0,p=1;
    for(int i=1;i<coeff.size();i++){
        res+=i*coeff[i]*p;
        p*=x;
    }
    return res;
}


pair<db,int> newtonRaphson(vector<db> &coeff, db x0, db tol, int maxIter){
    db x=x0;
    int iter;
    for(iter=1;iter<=maxIter;iter++){
        db fx=evalPoly(x,coeff);
        db fdx=evalPolyDerivative(x,coeff);
        if(fabs(fdx)<1e-12) break;
        db x1=x-fx/fdx;
        if(fabs(x1-x)<tol) return {x1,iter};
        x=x1;
    }
    return {x,iter};
}

int main(){
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int t; fin>>t;
    fout<<"Total Test Cases: "<<t<<"\n\n";

    for(int tc=1;tc<=t;tc++){
        int degree; fin>>degree;
        vector<db> coeff(degree+1);
        for(int i=0;i<=degree;i++) fin>>coeff[i];
        db start,end,tol; int maxIter;
        fin>>start>>end>>tol>>maxIter;

        fout<<"TEST CASE #"<<tc<<"\n";
        fout<<"Polynomial Degree: "<<degree<<"\n";
        fout<<"Given Polynomial: ";
        for(int i=degree;i>=0;i--){
            fout<<(int)coeff[i]<<"x^"<<i;
            if(i>0) fout<<" + ";
        }
        fout<<"\nInitial Guess Range: ["<<start<<", "<<end<<"]\n";
        fout<<"Error Tolerance: "<<tol<<"\n";
        fout<<"Max Iterations: "<<maxIter<<"\n";
        fout<<"Approximate Roots Found:\n";

        vector<pair<db,int>> roots;
        int steps=1000;
        db stepSize=(end-start)/steps;

        for(db x=start;x<end;x+=stepSize){
            db f1=evalPoly(x,coeff);
            db f2=evalPoly(x+stepSize,coeff);
            if(f1*f2<=0){
                db guess=(x+x+stepSize)/2;
                auto r=newtonRaphson(coeff,guess,tol,maxIter);
                bool unique=true;
                for(auto &p:roots){
                    if(fabs(p.first-r.first)<1e-6){unique=false; break;}
                }
                if(unique) roots.push_back(r);
            }
        }

        sort(roots.begin(),roots.end());
        for(int i=0;i<roots.size();i++){
            fout<<"  Root "<<i+1<<" = "<<fixed<<setprecision(6)<<roots[i].first
                <<" (Iterations: "<<roots[i].second<<")\n";
        }
        fout<<"\n";
    }

    fin.close(); fout.close();
    return 0;
}


```

#### Newton Raphson Input

```
3
2
1 -3 2
0 3 0.0001 100
3
1 0 -6 5
-2 4 0.0001 100
2
-2 5 -3
0 3 0.0001 100


```

#### Newton Raphson Output

```
Total Test Cases: 3

TEST CASE #1
Polynomial Degree: 2
Given Polynomial: 2x^2 + -3x^1 + 1x^0
Initial Guess Range: [0, 3]
Error Tolerance: 0.0001
Max Iterations: 100
Approximate Roots Found:
  Root 1 = 0.500000 (Iterations: 2)
  Root 2 = 1.000000 (Iterations: 2)

TEST CASE #2
Polynomial Degree: 3
Given Polynomial: 5x^3 + -6x^2 + 0x^1 + 1x^0
Initial Guess Range: [-2.000000, 4.000000]
Error Tolerance: 0.000100
Max Iterations: 100
Approximate Roots Found:
  Root 1 = -0.358258 (Iterations: 2)
  Root 2 = 0.558258 (Iterations: 2)
  Root 3 = 1.000000 (Iterations: 2)

TEST CASE #3
Polynomial Degree: 2
Given Polynomial: -3x^2 + 5x^1 + -2x^0
Initial Guess Range: [0.000000, 3.000000]
Error Tolerance: 0.000100
Max Iterations: 100
Approximate Roots Found:
  Root 1 = 0.666667 (Iterations: 2)
  Root 2 = 1.000000 (Iterations: 2)




```

---

### Secant Method

#### Secant Theory
 
The Secant Method is a numerical technique used to find the roots of a function. It is similar to the Newton-Raphson Method but does **not require the derivative** of the function. Instead, it uses a secant line through two initial guesses to approximate the root.
  
- Requires two initial guesses, `x0` and `x1`.  
- Does **not require the derivative** of the function.  
- Convergence is generally faster than Bisection but slower than Newton-Raphson.  
- Can fail if the secant line is nearly horizontal or if the function is highly nonlinear.

**Algorithm:**  
1. Choose two initial guesses `x0` and `x1`.  
2. Compute the next approximation using:  
x_{n+1} = x_n - f(x_n) * (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))
3. Replace `x_{n-1}` with `x_n` and `x_n` with `x_{n+1}`.  
4. Repeat step 2–3 until |f(x_n)| < tolerance or desired accuracy is achieved.

**Example:**  
Find a root of `f(x) = x^2 - 4` with initial guesses `x0 = 1`, `x1 = 3`:  
1. Compute next approximation:  
x2 = 3 - (3^2 - 4)(3 - 1)/((3^2 - 4) - (1^2 - 4)) = 3 - (52)/(5 - (-3)) = 3 - 10/8 = 1.75
2. Repeat:  
x3 = 1.75 - f(1.75)*(1.75 - 3)/(f(1.75) - f(3)) ≈ 2.0087
x4 ≈ 2.0000

Root converges to `x ≈ 2`.


#### Secant Code

```cpp
#include <bits/stdc++.h>
using namespace std;

typedef double db;


db evalPoly(vector<int>& coeff, db x) {
    db res = 0;
    int n = coeff.size();
    for(int i = 0; i < n; i++) {
        res += coeff[i] * pow(x, i);
    }
    return res;
}


pair<db,int> secantRoot(vector<int>& coeff, db x0, db x1, db tol, int maxIter) {
    db f0 = evalPoly(coeff, x0);
    db f1 = evalPoly(coeff, x1);
    db x2;
    int iter = 0;

    for(iter = 1; iter <= maxIter; iter++) {
        if(fabs(f1 - f0) < 1e-12) break; // avoid division by zero
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        if(fabs(evalPoly(coeff, x2)) < tol) return {x2, iter};
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = evalPoly(coeff, x1);
    }
    return {x2, iter};
}

int main() {
    ofstream fout("output.txt");
    fout << fixed << setprecision(6);

    int totalTestCases = 2;
    fout << "Total Test Cases: " << totalTestCases << "\n\n";


    vector<int> coeff1 = {1, -5, 6};
    db x0_1 = 0, x1_1 = 3, tol1 = 0.0001;
    int maxIter = 100;

    fout << "TEST CASE #1 :\n";
    fout << "Polynomial Degree: " << coeff1.size()-1 << "\n";
    fout << "Given Polynomial: ";
    for(int i = coeff1.size()-1; i>=0; i--) {
        fout << coeff1[i] << "x^" << i;
        if(i!=0) fout << " + ";
    }
    fout << "\nInitial Guesses: x0=" << x0_1 << ", x1=" << x1_1 << "\n";
    fout << "Root Bound: [" << x0_1 << ", " << x1_1 << "]\n";
    fout << "Error Tolerance: " << tol1 << "\n";
    fout << "Max Iterations: " << maxIter << "\n";


    vector<pair<db,int>> roots1;
    db step = (x1_1 - x0_1)/10.0;
    for(db a = x0_1; a < x1_1; a += step) {
        auto r = secantRoot(coeff1, a, a+step, tol1, maxIter);
        if(!isnan(r.first)) {
            bool unique = true;
            for(auto &p: roots1) if(fabs(p.first - r.first)<1e-4) unique=false;
            if(unique) roots1.push_back(r);
        }
    }

    for(int i=0;i<roots1.size();i++) {
        fout << "  Root " << i+1 << " in [" << roots1[i].first-step << ", " << roots1[i].first+step << "] = " << roots1[i].first << " (Iterations: " << roots1[i].second << ")\n";
    }
    fout << "\n";


    vector<int> coeff2 = {-3, 4, -1};
    db x0_2 = 0, x1_2 = 3, tol2 = 0.0001;

    fout << "TEST CASE #2 :\n";
    fout << "Polynomial Degree: " << coeff2.size()-1 << "\n";
    fout << "Given Polynomial: ";
    for(int i = coeff2.size()-1; i>=0; i--) {
        fout << coeff2[i] << "x^" << i;
        if(i!=0) fout << " + ";
    }
    fout << "\nInitial Guesses: x0=" << x0_2 << ", x1=" << x1_2 << "\n";
    fout << "Root Bound: [" << x0_2 << ", " << x1_2 << "]\n";
    fout << "Error Tolerance: " << tol2 << "\n";
    fout << "Max Iterations: " << maxIter << "\n";

    vector<pair<db,int>> roots2;
    step = (x1_2 - x0_2)/10.0;
    for(db a = x0_2; a < x1_2; a += step) {
        auto r = secantRoot(coeff2, a, a+step, tol2, maxIter);
        if(!isnan(r.first)) {
            bool unique = true;
            for(auto &p: roots2) if(fabs(p.first - r.first)<1e-4) unique=false;
            if(unique) roots2.push_back(r);
        }
    }

    for(int i=0;i<roots2.size();i++) {
        fout << "  Root " << i+1 << " in [" << roots2[i].first-step << ", " << roots2[i].first+step << "] = " << roots2[i].first << " (Iterations: " << roots2[i].second << ")\n";
    }

    fout.close();
    return 0;
}

```

#### Secant Input

```
2
6 -5 1
0 3
-3 4 -1
0 3

```

#### Secant Output

```
Total Test Cases: 2

TEST CASE #1 :
Polynomial Degree: 2
Given Polynomial: 6x^2 + -5x^1 + 1x^0
Initial Guesses: x0=0.000000, x1=3.000000
Root Bound: [0.000000, 3.000000]
Error Tolerance: 0.000100
Max Iterations: 100
  Root 1 in [0.033327, 0.633327] = 0.333327 (Iterations: 4)
  Root 2 in [0.200076, 0.800076] = 0.500076 (Iterations: 5)

TEST CASE #2 :
Polynomial Degree: 2
Given Polynomial: -1x^2 + 4x^1 + -3x^0
Initial Guesses: x0=0.000000, x1=3.000000
Root Bound: [0.000000, 3.000000]
Error Tolerance: 0.000100
Max Iterations: 100
  Root 1 in [0.700000, 1.300000] = 1.000000 (Iterations: 5)
  Root 2 in [2.699993, 3.299993] = 2.999993 (Iterations: 6)

```

---

## Differential Equation Solving

---

### Runge-Kutta 4th Order Method

#### Runge-Kutta Theory

 
The Runge-Kutta 4th Order Method (RK4) is a numerical technique used to solve ordinary differential equations (ODEs) of the form `dy/dx = f(x, y)` with a given initial condition `y(x0) = y0`. RK4 is widely used because it provides high accuracy without requiring higher derivatives of the function.

- Provides **fourth-order accuracy**, meaning the error per step is proportional to `h^5` and cumulative error is proportional to `h^4`.  
- Iterative method that uses **four weighted slopes** to estimate the next value.  
- More accurate than Euler’s method for the same step size.  
- Requires computation of `f(x, y)` four times per step.

**Algorithm:**  
Given `dy/dx = f(x, y)` and step size `h`:

k1 = h * f(x_n, y_n)
k2 = h * f(x_n + h/2, y_n + k1/2)
k3 = h * f(x_n + h/2, y_n + k2/2)
k4 = h * f(x_n + h, y_n + k3)
y_{n+1} = y_n + (1/6)(k1 + 2k2 + 2*k3 + k4)
x_{n+1} = x_n + h


**Example:**  
Solve `dy/dx = x + y`, with `y(0) = 1` and step size `h = 0.1`:  
1. Compute slopes:  
k1 = 0.1 * f(0, 1) = 0.1*(0+1) = 0.1
k2 = 0.1 * f(0.05, 1+0.05) = 0.1*(0.05+1.05) = 0.11
k3 = 0.1 * f(0.05, 1+0.055) = 0.1*(0.05+1.055) = 0.1105
k4 = 0.1 * f(0.1, 1+0.1105) = 0.1*(0.1+1.1105) = 0.12105
2. Compute next value:  
y1 = 1 + (1/6)(0.1 + 20.11 + 2*0.1105 + 0.12105) ≈ 1.11034




#### Runge-Kutta Code

```cpp
#include<bits/stdc++.h>
using namespace std;

double f(double x,double y){
    return (x+y)/2.0;
}

int main(){
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int t;
    fin>>t;

    for(int tc=1;tc<=t;tc++){
        double x0,y0,xn,h;
        fin>>x0>>y0>>xn>>h;

        int n=(xn-x0)/h;
        double x=x0;
        double y=y0;

        for(int i=1;i<=n;i++){
            double k1=h*f(x,y);
            double k2=h*f(x+h/2.0,y+k1/2.0);
            double k3=h*f(x+h/2.0,y+k2/2.0);
            double k4=h*f(x+h,y+k3);

            y=y+(k1+2*k2+2*k3+k4)/6.0;
            x=x+h;
        }

        fout<<fixed<<setprecision(3);
        fout<<"TEST CASE #"<<tc<<"\n\n";
        fout<<"Initial x0: "<<x0<<", y0: "<<y0<<"\n";
        fout<<"Final x: "<<xn<<"\n";
        fout<<"Step size (h): "<<h<<"\n";
        fout<<"Number of steps: "<<n<<"\n";
        fout<<setprecision(6);
        fout<<"Result: y("<<xn<<") = "<<y<<"\n\n";
    }

    fin.close();
    fout.close();
    return 0;
}

```

#### Runge-Kutta Input

```
3
0 1 1.5 0.01
0 2 2.0 0.001
1 1 2.5 0.005

```

#### Runge-Kutta Output

```
TEST CASE #1

Initial x0: 0.000, y0: 1.000
Final x: 1.500
Step size (h): 0.010
Number of steps: 150
Result: y(1.500000) = 2.851000

TEST CASE #2

Initial x0: 0.000, y0: 2.000
Final x: 2.000
Step size (h): 0.001
Number of steps: 2000
Result: y(2.000000) = 6.873127

TEST CASE #3

Initial x0: 1.000, y0: 1.000
Final x: 2.500
Step size (h): 0.005
Number of steps: 300
Result: y(2.500000) = 3.968000


```

---

## Interpolation Methods

---

### Newton Forward Interpolation

#### Newton Forward Theory

 
Newton Forward Interpolation is used to estimate the value of a function at a point near the **beginning** of a table of equally spaced data points.
  
- Works with **equally spaced data** points.  
- Uses the **forward difference operator (Δ)**.  
- More accurate when the value to estimate is **near the first data point**.  
- Iterative polynomial construction using forward differences.

**Algorithm:**  
f(x) ≈ f(x0) + uΔf(x0) + u(u-1)/2! Δ²f(x0) + ... + u(u-1)...(u-n+1)/n! Δⁿf(x0)

Where:  
- `u = (x - x0)/h`, `h` is the spacing between x-values  
- Δ is the forward difference operator  
- `n` is the number of terms used

**Example:**  
Given data: x = 0,1,2; f(x) = 1,2,5  
Estimate f(0.5):  
u = (0.5 - 0)/1 = 0.5
f(0.5) ≈ f(0) + 0.5Δf(0) + 0.5(0.5-1)/2 * Δ²f(0)



#### Newton Forward Code

```cpp
#include <bits/stdc++.h>
using namespace std;

double factorial(int n){
    double f=1;
    for(int i=2;i<=n;i++) f*=i;
    return f;
}

vector<vector<double>> forwardTable(const vector<double>& y){
    int n=y.size();
    vector<vector<double>> table(n, vector<double>(n,0));
    for(int i=0;i<n;i++) table[i][0]=y[i];
    for(int j=1;j<n;j++)
        for(int i=0;i<n-j;i++)
            table[i][j]=table[i+1][j-1]-table[i][j-1];
    return table;
}

double forwardInterpolation(const vector<double>& x, const vector<vector<double>>& t, double p){
    double h=x[1]-x[0], u=(p-x[0])/h, res=t[0][0], term=1;
    for(int i=1;i<t.size();i++){
        term*=(u-i+1);
        res+=term*t[0][i]/factorial(i);
    }
    return res;
}

int main(){
    ifstream in("input.txt");
    ofstream out("output.txt");
    int tests; in>>tests;
    out<<fixed<<setprecision(6);
    for(int tc=1;tc<=tests;tc++){
        int n; in>>n;
        vector<double> x(n), y(n);
        for(int i=0;i<n;i++) in>>x[i];
        for(int i=0;i<n;i++) in>>y[i];
        double p; in>>p;
        auto t=forwardTable(y);
        out<<"Test Case "<<tc<<"\n";
        out<<"Data Points:\n";
        for(int i=0;i<n;i++) out<<"("<<x[i]<<", "<<y[i]<<")\n";
        out<<"Interpolation Point: "<<p<<"\nStep Size: "<<x[1]-x[0]<<"\n";
        out<<"Forward Difference Table:\n";
        for(int j=0;j<n;j++){
            for(int i=0;i<n-j;i++) out<<setw(6)<<t[i][j]<<" ";
            out<<"\n";
        }
        out<<"Interpolated Value: "<<forwardInterpolation(x,t,p)<<"\n\n";
    }
}

```

#### Newton Forward Input

```
2
5
0 1 2 3 4
1 2 4 8 16
2.5
4
1 2 3 4
1 8 27 64
2.5

```

#### Newton Forward Output

```
Test Case 1
Data Points:
(0, 1)
(1, 2)
(2, 4)
(3, 8)
(4, 16)
Interpolation Point: 2.5
Step Size: 1
Forward Difference Table:
     1      1      1      1      1 
     2      2      2      2 
     4      4      4 
     8      8 
    16 
Interpolated Value: 5.648438

Test Case 2
Data Points:
(1, 1)
(2, 8)
(3, 27)
(4, 64)
Interpolation Point: 2.5
Step Size: 1
Forward Difference Table:
     1      7     12      6 
     8     19     18 
    27     37 
    64 
Interpolated Value: 15.625000

```

---

### Newton Backward Interpolation

#### Newton Backward Theory

 
Newton Backward Interpolation is used to estimate the value of a function at a point near the **end** of a table of equally spaced data points.
 
- Works with **equally spaced data** points.  
- Uses the **backward difference operator (∇)**.  
- More accurate when the value to estimate is **near the last data point**.  
- Iterative polynomial construction using backward differences.

**Algorithm:**  
f(x) ≈ f(xn) + u∇f(xn) + u(u+1)/2! ∇²f(xn) + ... + u(u+1)...(u+n-1)/n! ∇ⁿf(xn)

Where:  
- `u = (x - xn)/h`, `h` is the spacing between x-values  
- ∇ is the backward difference operator  
- `n` is the number of terms used

**Example:**  
Given data: x = 0,1,2; f(x) = 1,2,5  
Estimate f(1.5):  
u = (1.5 - 2)/1 = -0.5
f(1.5) ≈ f(2) + (-0.5)*∇f(2) + (-0.5)(-0.5+1)/2 * ∇²f(2)



#### Newton Backward Code

```cpp
#include <bits/stdc++.h>
using namespace std;

double factorial(int n){
    double f=1;
    for(int i=2;i<=n;i++) f*=i;
    return f;
}

vector<vector<double>> backwardTable(const vector<double>& y){
    int n=y.size();
    vector<vector<double>> t(n, vector<double>(n,0));
    for(int i=0;i<n;i++) t[i][0]=y[i];
    for(int j=1;j<n;j++){
        for(int i=n-1;i>=j;i--) t[i][j]=t[i][j-1]-t[i-1][j-1];
    }
    return t;
}

double backwardInterpolation(const vector<double>& x,const vector<vector<double>>& t,double p){
    int n=x.size();
    double h=x[1]-x[0], u=(p-x[n-1])/h, res=t[n-1][0], term=1;
    for(int i=1;i<n;i++){
        term*=(u+i-1);
        res+=term*t[n-1][i]/factorial(i);
    }
    return res;
}

int main(){
    ifstream in("input.txt");
    ofstream out("output.txt");
    int tests; in>>tests;
    out<<fixed<<setprecision(6);
    for(int tc=1;tc<=tests;tc++){
        int n; in>>n;
        vector<double> x(n), y(n);
        for(int i=0;i<n;i++) in>>x[i];
        for(int i=0;i<n;i++) in>>y[i];
        double p; in>>p;
        auto t=backwardTable(y);
        out<<"Test Case "<<tc<<"\n";
        out<<"Data Points:\n";
        for(int i=0;i<n;i++) out<<"("<<x[i]<<", "<<y[i]<<")\n";
        out<<"Interpolation Point: "<<p<<"\nStep Size: "<<x[1]-x[0]<<"\n";
        out<<"Backward Difference Table:\n";
        for(int j=0;j<n;j++){
            for(int i=j;i<n;i++) out<<setw(6)<<t[i][j]<<" ";
            out<<"\n";
        }
        out<<"Interpolated Value: "<<backwardInterpolation(x,t,p)<<"\n\n";
    }
}

```

#### Newton Backward Input

```
2
5
0 1 2 3 4
1 2 4 8 16
2.5
4
1 2 3 4
1 8 27 64
2.5

```

#### Newton Backward Output

```
Test Case 1
Data Points:
(0, 1)
(1, 2)
(2, 4)
(3, 8)
(4, 16)
Interpolation Point: 2.5
Step Size: 1
Backward Difference Table:
     1      1      1      1      1 
     2      2      2      2 
     4      4      4 
     8      8 
    16 
Interpolated Value: 5.648438

Test Case 2
Data Points:
(1, 1)
(2, 8)
(3, 27)
(4, 64)
Interpolation Point: 2.5
Step Size: 1
Backward Difference Table:
     1      7     12      6 
     8     19     18 
    27     37 
    64 
Interpolated Value: 15.625000

```

---

### Newton Divided Difference Interpolation

#### Newton Divided Difference Theory

  
Newton Divided Difference Interpolation is used for **unequally spaced data points**. It constructs an interpolation polynomial using divided differences.

- Works for **unequally spaced data**.  
- Uses **divided differences** instead of forward/backward differences.  
- Flexible for any point in the dataset.  
- Polynomial construction avoids the limitation of equally spaced data.

**Algorithm:**  
f(x) = f[x0] + (x - x0)f[x0,x1] + (x - x0)(x - x1)f[x0,x1,x2] + ... + (x - x0)...(x - xn-1)f[x0,x1,...,xn]

Where:  
- `f[x0,x1,...,xn]` are the divided differences  
- Works for both equally and unequally spaced data

**Example:**  
Given data: x = 1,2,4; f(x) = 1,3,15  
Estimate f(3):  
Compute divided differences f[x0,x1], f[x1,x2], f[x0,x1,x2], ...
f(3) ≈ f[x0] + (3 - 1)f[x0,x1] + (3 - 1)(3 - 2)*f[x0,x1,x2]



#### Newton Divided Difference Code

```cpp
#include <bits/stdc++.h>
using namespace std;

pair<vector<vector<double>>, int> dividedTable(const vector<double>& x, vector<double> y) {
    int n = y.size();
    vector<vector<double>> table(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) table[i][0] = y[i];
    for (int j = 1; j < n; j++)
        for (int i = 0; i < n - j; i++)
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x[i + j] - x[i]);
    int order = 0;
    for (int i = 0; i < n; i++)
        if (abs(table[0][i]) > 1e-6) order = i;
    return {table, order};
}

double interpolate(const vector<double>& x, const vector<vector<double>>& table, double p) {
    double res = table[0][0], term = 1;
    for (int i = 1; i < table.size(); i++) {
        term *= (p - x[i - 1]);
        res += table[0][i] * term;
    }
    return res;
}

int main() {
    ifstream in("input.txt");
    ofstream out("output.txt");
    int tests;
    in >> tests;
    out << fixed << setprecision(6);
    for (int tc = 1; tc <= tests; tc++) {
        int n; in >> n;
        vector<double> x(n), y(n);
        for (int i = 0; i < n; i++) in >> x[i];
        for (int i = 0; i < n; i++) in >> y[i];
        double p; in >> p;
        auto [table, order] = dividedTable(x, y);
        double val = interpolate(x, table, p);
        double error = (table[0][order + 1] != 0 ? table[0][order + 1] : 0);
        out << "Test Case " << tc << "\n";
        out << "Data Points:\n";
        for (int i = 0; i < n; i++)
            out << "(" << x[i] << ", " << y[i] << ")\n";
        out << "Interpolation Point: " << p << "\n";
        out << "Divided Difference Table:\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n - i; j++)
                out << table[i][j] << " ";
            out << "\n";
        }
        out << "Polynomial Order Detected: " << order << "\n";
        out << "Interpolated Value: " << val << "\n";
        out << "Estimated Error (next term): " << error << "\n\n";
    }
}

```

#### Newton Divided Difference Input

```
2
5
0 1 2 3 4
1 2 4 8 16
2.5
4
1 2 3 4
1 8 27 64
2.5

```

#### Newton Divided Difference Output

```
Test Case 1
Data Points:
(0, 1)
(1, 2)
(2, 4)
(3, 8)
(4, 16)
Interpolation Point: 2.5
Divided Difference Table:
1 1 1 0 0 
2 1 2 4 
4 4 8 
8 12 
16 
Polynomial Order Detected: 4
Interpolated Value: 5.648438
Estimated Error (next term): 0

Test Case 2
Data Points:
(1, 1)
(2, 8)
(3, 27)
(4, 64)
Interpolation Point: 2.5
Divided Difference Table:
1 7 6 1 
8 19 18 
27 37 
64 
Polynomial Order Detected: 3
Interpolated Value: 15.625000
Estimated Error (next term): 0

```

---

## Numerical Differentiation

---

### Differentiation by Forward Interpolation

#### Differentiation Forward Theory

### Algorithm
1. Choose equally spaced data near the start.
2. Construct the forward difference table.
3. Compute step size h.
4. Find u = (x − x0) / h.
5. Apply forward differentiation formula.

### Mathematical Formula

First derivative:

f'(x) = [ Δf0 + (2u − 1)Δ²f0/2 + (3u² − 6u + 2)Δ³f0/6 + … ] / h

Second derivative:

f''(x) = [ Δ²f0 + (u − 1)Δ³f0 + … ] / h²

where:

u = (x − x0) / h

### Brief
Forward differentiation is used near the beginning of data.
It requires equally spaced points.
The method is simple and systematic.
Accuracy improves with higher differences.
It is suitable for polynomial functions.


#### Differentiation Forward Code

```cpp
#include<bits/stdc++.h>
using namespace std;

double f(double x)
{
    return x*x*x;
}

void forwardDiff(ifstream &in,ofstream &out)
{
    int T;
    in>>T;

    for(int c=1;c<=T;c++)
    {
        int i;
        double a,b,p;
        in>>i;
        in>>a>>b;
        in>>p;

        int n=i+1;
        double h=(b-a)/i;

        vector<double> x(n),y(n);
        for(int k=0;k<n;k++)
        {
            x[k]=a+k*h;
            y[k]=f(x[k]);
        }

        vector<vector<double>> d(n,vector<double>(n));
        for(int k=0;k<n;k++)
            d[k][0]=y[k];

        for(int j=1;j<n;j++)
            for(int k=0;k<n-j;k++)
                d[k][j]=d[k+1][j-1]-d[k][j-1];

        double u=(p-a)/h;

        double f1=d[0][1];
        if(n>2) f1+=(2*u-1)*d[0][2]/2.0;
        if(n>3) f1+=(3*u*u-6*u+2)*d[0][3]/6.0;
        f1/=h;

        double f2=0;
        if(n>2) f2+=d[0][2];
        if(n>3) f2+=(u-1)*d[0][3];
        f2/=(h*h);

        out<<"Case "<<c<<":\n";
        out<<"Forward Difference Table:\n";
        for(int k=0;k<n;k++)
        {
            for(int j=0;j<n-k;j++)
                out<<fixed<<setprecision(3)<<d[k][j]<<" ";
            out<<"\n";
        }

        out<<"First Derivative at p = "<<fixed<<setprecision(3)<<p<<" : "<<fixed<<setprecision(3)<<f1<<"\n";
        out<<"Second Derivative at p = "<<fixed<<setprecision(3)<<p<<" : "<<fixed<<setprecision(3)<<f2<<"\n\n";
    }
}

int main()
{
    ifstream in("diff_input.txt");
    ofstream out("diff_output.txt");
    forwardDiff(in,out);
}

```

#### Differentiation Forward Input

```
3
4
0 4
1.5
4
0 4
2.5
5
0 5
1.2

```

#### Differentiation Forward Output

```
Case 1:
Forward Difference Table:
0.000 1.000 6.000 6.000 0.000
1.000 7.000 12.000 6.000
8.000 19.000 18.000
27.000 37.000
64.000
First Derivative at p = 1.500 : 6.750
Second Derivative at p = 1.500 : 9.000

Case 2:
Forward Difference Table:
0.000 1.000 6.000 6.000 0.000
1.000 7.000 12.000 6.000
8.000 19.000 18.000
27.000 37.000
64.000
First Derivative at p = 2.500 : 18.750
Second Derivative at p = 2.500 : 15.000

Case 3:
Forward Difference Table:
0.000 1.000 6.000 6.000 0.000 0.000
1.000 7.000 12.000 6.000 0.000
8.000 19.000 18.000 6.000
27.000 37.000 24.000
64.000 61.000
125.000
First Derivative at p = 1.200 : 4.320
Second Derivative at p = 1.200 : 7.200

```

---

### Differentiation by Backward Interpolation

#### Differentiation Backward Theory

### Algorithm
1. Choose equally spaced data near the end.
2. Construct the backward difference table.
3. Compute step size h.
4. Find s = (x − xn) / h.
5. Apply backward differentiation formula.

### Mathematical Formula

First derivative:

f'(x) = [ ∇f_n + (2s + 1)∇²f_n/2 + (3s² + 6s + 2)∇³f_n/6 + … ] / h

where:

s = (x − x_n) / h

### Brief
Backward differentiation is used near the last data point.
It uses backward differences.
The method gives good accuracy near the end.
It is effective for smooth functions.
Backward interpolation complements forward interpolation.

#### Differentiation Backward Code

```cpp
#include<bits/stdc++.h>
using namespace std;

double f(double x,vector<double> c)
{
    double r=0,p=1;
    for(double v:c)
    {
        r+=v*p;
        p*=x;
    }
    return r;
}

double df(double x,vector<double> c)
{
    double r=0,p=1;
    for(int i=1;i<(int)c.size();i++)
    {
        r+=i*c[i]*p;
        p*=x;
    }
    return r;
}

void backwardDiff(ifstream &in,ofstream &out)
{
    int T;
    in>>T;

    for(int tc=1;tc<=T;tc++)
    {
        int n;
        in>>n;

        vector<double> c(n+1);
        for(int i=0;i<=n;i++) in>>c[i];

        int p;
        in>>p;

        vector<double> x(p),y(p);
        for(int i=0;i<p;i++)
        {
            in>>x[i];
            y[i]=f(x[i],c);
        }

        double d;
        in>>d;

        double h=x[1]-x[0];

        vector<vector<double>> bd(p,vector<double>(p));
        for(int i=0;i<p;i++) bd[i][0]=y[i];

        for(int j=1;j<p;j++)
            for(int i=p-1;i>=j;i--)
                bd[i][j]=bd[i][j-1]-bd[i-1][j-1];

        double s=(d-x[p-1])/h;

        double num=bd[p-1][1];
        if(p>2) num+=(2*s+1)*bd[p-1][2]/2.0;
        if(p>3) num+=(3*s*s+6*s+2)*bd[p-1][3]/6.0;
        num/=h;

        double exact=df(d,c);
        double err=fabs((exact-num)/exact)*100.0;

        out<<"Case "<<tc<<":\n";
        out<<"Backward Difference Table:\n";
        for(int i=0;i<p;i++)
        {
            for(int j=0;j<=i;j++)
                out<<fixed<<setprecision(3)<<bd[i][j]<<" ";
            out<<"\n";
        }

        out<<"Numerical Derivative : "<<fixed<<setprecision(3)<<num<<"\n";
        out<<"Exact Derivative     : "<<fixed<<setprecision(3)<<exact<<"\n";
        out<<"Error Percentage     : "<<fixed<<setprecision(3)<<err<<" %\n\n";
    }
}

int main()
{
    ifstream in("backward_input.txt");
    ofstream out("backward_output.txt");
    backwardDiff(in,out);
}

```

#### Differentiation Backward Input

```
3
3
0 0 0 1
5
0 1 2 3 4
3.8
2
1 0 1
4
1 2 3 4
3.6
3
1 -2 0 1
5
0 1 2 3 4
3.5

```

#### Differentiation Backward Output

```
Case 1:
Backward Difference Table:
0.000
1.000 1.000
8.000 7.000 6.000
27.000 19.000 12.000 6.000
64.000 37.000 18.000 6.000 0.000
Numerical Derivative : 43.320
Exact Derivative     : 43.320
Error Percentage     : 0.000 %

Case 2:
Backward Difference Table:
2.000
5.000 3.000
10.000 5.000 2.000
17.000 7.000 2.000 0.000
Numerical Derivative : 7.200
Exact Derivative     : 7.200
Error Percentage     : 0.000 %

Case 3:
Backward Difference Table:
1.000
0.000 -1.000
5.000 5.000 6.000
22.000 17.000 12.000 6.000
53.000 31.000 14.000 2.000 -4.000
Numerical Derivative : 33.750
Exact Derivative     : 34.750
Error Percentage     : 2.877 %

```

---

## Curve Fitting / Regression

---

### Linear Regression

#### Linear Regression Theory

  
Linear Regression is a statistical method used to model the relationship between a dependent variable `y` and an independent variable `x` by fitting a straight line (`y = mx + c`) to the data points.
  
- Models the data with a **straight line**.  
- Parameters `m` (slope) and `c` (intercept) are determined using the **least squares method**.  
- Minimizes the sum of the squares of the vertical deviations of the points from the line.  
- Suitable when the relationship between variables is approximately linear.

**Algorithm:**  
m = [nΣ(xy) - ΣxΣy] / [nΣ(x^2) - (Σx)^2]
c = [Σy - mΣx] / n

Where:  
- `n` = number of data points  
- Σx, Σy = sum of x and y values  
- Σ(xy) = sum of products of x and y

**Example:**  
Given data points: (1,2), (2,3), (3,5)  
m = [3*(12 + 23 + 35) - (1+2+3)(2+3+5)] / [3*(1^2+2^2+3^2) - (1+2+3)^2] = 1.5
c = [2+3+5 - 1.5*(1+2+3)] / 3 = 0.33
y ≈ 1.5x + 0.33



#### Linear Regression Code

```cpp
#include <bits/stdc++.h>
using namespace std;

pair<double,double> linearRegression(const vector<double>& x,const vector<double>& y){
    int n=x.size();
    double sx=0,sy=0,sxy=0,sx2=0;
    for(int i=0;i<n;i++){
        sx+=x[i];
        sy+=y[i];
        sxy+=x[i]*y[i];
        sx2+=x[i]*x[i];
    }
    double b=(n*sxy-sx*sy)/(n*sx2-sx*sx);
    double a=(sy-b*sx)/n;
    return {a,b};
}

int main(){
    ifstream in("input.txt");
    ofstream out("output.txt");

    int tests;
    in>>tests;
    out<<fixed<<setprecision(6);

    for(int t=1;t<=tests;t++){
        int n;
        in>>n;
        vector<double> x(n),y(n);
        for(int i=0;i<n;i++) in>>x[i];
        for(int i=0;i<n;i++) in>>y[i];

        auto result=linearRegression(x,y);

        out<<"Linear Regression\n";
        out<<"Test Case "<<t<<"\n";
        out<<"Number of Points: "<<n<<"\n";
        out<<"Data Points:\n";
        for(int i=0;i<n;i++) out<<"("<<x[i]<<", "<<y[i]<<")\n";
        out<<"Target Equation: y = a + bx\n";
        out<<"Intercept (a): "<<result.first<<"\n";
        out<<"Slope (b): "<<result.second<<"\n";
        out<<"Linear Relationship Equation: y = "<<result.first<<" + "<<result.second<<"x\n\n";
    }
}

```

#### Linear Regression Input

```
3
7
1 2 3 4 5 6 7
3 4 4 5 8 9 10
5
1 2 3 4 5
2 4 6 8 10
6
1 2 3 4 5 6
5 7 9 11 13 15

```

#### Linear Regression Output

```
Linear Regression
Test Case 1
Number of Points: 7
Data Points:
(1, 3)
(2, 4)
(3, 4)
(4, 5)
(5, 8)
(6, 9)
(7, 10)
Target Equation: y = a + bx
Intercept (a): 1.095238
Slope (b): 1.190476
Linear Relationship Equation: y = 1.095238 + 1.190476x

Linear Regression
Test Case 2
Number of Points: 5
Data Points:
(1, 2)
(2, 4)
(3, 6)
(4, 8)
(5, 10)
Target Equation: y = a + bx
Intercept (a): 0.000000
Slope (b): 2.000000
Linear Relationship Equation: y = 0 + 2x

Linear Regression
Test Case 3
Number of Points: 6
Data Points:
(1, 5)
(2, 7)
(3, 9)
(4, 11)
(5, 13)
(6, 15)
Target Equation: y = a + bx
Intercept (a): 3.000000
Slope (b): 2.000000
Linear Relationship Equation: y = 3 + 2x

```

---

### Polynomial Regression

#### Polynomial Regression Theory

  
Polynomial Regression is a regression technique where the relationship between the independent variable `x` and the dependent variable `y` is modeled as an **nth-degree polynomial**.

- Fits data with a **curve** instead of a straight line.  
- Useful when data shows a **non-linear trend**.  
- Parameters are estimated using **least squares method**.  
- Degree of the polynomial is chosen based on data complexity.

**Algorithm:**  
y = a0 + a1x + a2x^2 + ... + an*x^n

Where:  
- `n` = degree of the polynomial  
- `a0, a1, ..., an` = coefficients determined using least squares

**Example:**  
Given data: x = 1,2,3; y = 1,4,9  
Fit a quadratic polynomial: y ≈ a0 + a1x + a2x^2
Solve simultaneous equations to find coefficients: a0=0, a1=0, a2=1
y ≈ x^2



#### Polynomial Regression Code

```cpp
#include <bits/stdc++.h>
using namespace std;

vector<double> gaussianElimination(vector<vector<double>> A, vector<double> B){
    int n=B.size();
    for(int i=0;i<n;i++){
        for(int k=i+1;k<n;k++){
            double factor=A[k][i]/A[i][i];
            for(int j=i;j<n;j++) A[k][j]-=factor*A[i][j];
            B[k]-=factor*B[i];
        }
    }
    vector<double> X(n);
    for(int i=n-1;i>=0;i--){
        X[i]=B[i];
        for(int j=i+1;j<n;j++) X[i]-=A[i][j]*X[j];
        X[i]/=A[i][i];
    }
    return X;
}

vector<double> polynomialRegression(const vector<double>& x,const vector<double>& y,int degree){
    int n=degree+1,m=x.size();
    vector<vector<double>> A(n,vector<double>(n,0));
    vector<double> B(n,0);

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            for(int k=0;k<m;k++)
                A[i][j]+=pow(x[k],i+j);
        for(int k=0;k<m;k++)
            B[i]+=y[k]*pow(x[k],i);
    }
    return gaussianElimination(A,B);
}

int main(){
    ifstream in("input.txt");
    ofstream out("output.txt");
    out<<fixed<<setprecision(6);

    int tests;
    in>>tests;

    for(int t=1;t<=tests;t++){
        int n,degree;
        in>>n>>degree;
        vector<double> x(n),y(n);
        for(int i=0;i<n;i++) in>>x[i];
        for(int i=0;i<n;i++) in>>y[i];

        auto c=polynomialRegression(x,y,degree);

        out<<"Polynomial Regression\n";
        out<<"Test Case "<<t<<"\n";
        out<<"Polynomial Degree: "<<degree<<"\n";
        out<<"Target Equation: y = a0 + a1x + a2x^2 + ... + a"<<degree<<"x^"<<degree<<"\n";
        out<<"Data Points:\n";
        for(int i=0;i<n;i++) out<<"("<<x[i]<<", "<<y[i]<<")\n";
        out<<"Coefficients:\n";
        for(int i=0;i<c.size();i++) out<<"a"<<i<<" = "<<c[i]<<"\n";
        out<<"Curve Equation: y = ";
        for(int i=0;i<c.size();i++){
            if(i>0) out<<" + ";
            out<<c[i];
            if(i>0) out<<"x^"<<i;
        }
        out<<"\n\n";
    }
}

```

#### Polynomial Regression Input

```
3
5 2
1 2 3 4 5
6 11 18 27 38
6 3
1 2 3 4 5 6
2 9 28 65 126 217
5 1
1 2 3 4 5
3 5 7 9 11

```

#### Polynomial Regression Output

```
Polynomial Regression
Test Case 1
Polynomial Degree: 2
Target Equation: y = a0 + a1x + a2x^2 + ... + a2x^2
Data Points:
(1, 6)
(2, 11)
(3, 18)
(4, 27)
(5, 38)
Coefficients:
a0 = 1.000000
a1 = 2.000000
a2 = 1.000000
Curve Equation: y = 1 + 2x^1 + 1x^2

Polynomial Regression
Test Case 2
Polynomial Degree: 3
Target Equation: y = a0 + a1x + a2x^2 + ... + a3x^3
Data Points:
(1, 2)
(2, 9)
(3, 28)
(4, 65)
(5, 126)
(6, 217)
Coefficients:
a0 = 0.000000
a1 = 0.000000
a2 = 0.000000
a3 = 1.000000
Curve Equation: y = 1x^3

Polynomial Regression
Test Case 3
Polynomial Degree: 1
Target Equation: y = a0 + a1x + ... + a1x^1
Data Points:
(1, 3)
(2, 5)
(3, 7)
(4, 9)
(5, 11)
Coefficients:
a0 = 1.000000
a1 = 2.000000
Curve Equation: y = 1 + 2x

```

---

### Transcendental Regression

#### Transcendental Regression Theory


 
Transcendental Regression is used to fit data using **non-polynomial functions** such as exponential, logarithmic, or trigonometric functions. It is suitable when the relationship between variables is not polynomial.

- Fits data to functions like **y = a*e^(bx)**, **y = a*ln(x) + b**, or **y = a*sin(bx) + c**.  
- Parameters are determined using **least squares or nonlinear optimization**.  
- Useful for data with **exponential growth, decay, or periodic patterns**.

**Algorithm:**  
Example: For exponential regression:  
y = a * e^(bx)
Take ln(y) = ln(a) + bx => Apply linear regression to ln(y) vs x


**Example:**  
Given data: x = 0,1,2; y = 2,4.4,9.6  
Assume y = ae^(bx)
ln(y) = ln(a) + b*x
Apply linear regression on ln(y) vs x to find ln(a) and b
Compute a = e^(ln(a)), b



#### Transcendental Regression Code

```cpp
#include <bits/stdc++.h>
using namespace std;

pair<double,double> linearFit(const vector<double>& X,const vector<double>& Y){
    int n=X.size();
    double sx=0,sy=0,sxy=0,sx2=0;
    for(int i=0;i<n;i++){
        sx+=X[i];
        sy+=Y[i];
        sxy+=X[i]*Y[i];
        sx2+=X[i]*X[i];
    }
    double b=(n*sxy-sx*sy)/(n*sx2-sx*sx);
    double a=(sy-b*sx)/n;
    return {a,b};
}

int main(){
    ifstream in("input.txt");
    ofstream out("output.txt");
    out<<fixed<<setprecision(6);

    int tests;
    in>>tests;

    for(int t=1;t<=tests;t++){
        int model,n;
        in>>model>>n;

        vector<double> x(n),y(n);
        for(int i=0;i<n;i++) in>>x[i];
        for(int i=0;i<n;i++) in>>y[i];

        out<<"Transcendental Regression\n";
        out<<"Test Case "<<t<<"\n";
        out<<"Data Points:\n";
        for(int i=0;i<n;i++)
            out<<"("<<x[i]<<", "<<y[i]<<")\n";

        if(model==1){
            vector<double> X(n),Y(n);
            for(int i=0;i<n;i++){
                X[i]=x[i];
                Y[i]=log(y[i]);
            }
            auto r=linearFit(X,Y);
            out<<"Model: y = a e^(bx)\n";
            out<<"a: "<<exp(r.first)<<"\n";
            out<<"b: "<<r.second<<"\n";
            out<<"Curve Equation: y = "<<exp(r.first)<<" e^("<<r.second<<"x)\n";
        }

        if(model==2){
            vector<double> X(n),Y(n);
            for(int i=0;i<n;i++){
                X[i]=log(x[i]);
                Y[i]=log(y[i]);
            }
            auto r=linearFit(X,Y);
            out<<"Model: y = a x^b\n";
            out<<"a: "<<exp(r.first)<<"\n";
            out<<"b: "<<r.second<<"\n";
            out<<"Curve Equation: y = "<<exp(r.first)<<" x^"<<r.second<<"\n";
        }

        if(model==3){
            vector<double> X(n),Y(n);
            for(int i=0;i<n;i++){
                X[i]=exp(x[i]/4);
                Y[i]=y[i];
            }
            auto r=linearFit(X,Y);
            out<<"Model: y = a + b e^(x/4)\n";
            out<<"a: "<<r.first<<"\n";
            out<<"b: "<<r.second<<"\n";
            out<<"Curve Equation: y = "<<r.first<<" + "<<r.second<<" e^(x/4)\n";
        }

        out<<"\n";
    }
}

```

#### Transcendental Regression Input

```
3
1 5
1 2 3 4 5
3 8 20 55 148
2 5
1 2 3 4 5
2 8 18 32 50
3 5
1 2 3 4 5
6 9 14 22 35

```

#### Transcendental Regression Output

```
Transcendental Regression
Test Case 1
Data Points:
(1, 3)
(2, 8)
(3, 20)
(4, 55)
(5, 148)
Model: y = a e^(bx)
a: 1.523602
b: 0.998317
Curve Equation: y = 1.523602 e^(0.998317x)

Transcendental Regression
Test Case 2
Data Points:
(1, 2)
(2, 8)
(3, 18)
(4, 32)
(5, 50)
Model: y = a x^b
a: 1.997814
b: 2.000128
Curve Equation: y = 1.997814 x^2.000128

Transcendental Regression
Test Case 3
Data Points:
(1, 6)
(2, 9)
(3, 14)
(4, 22)
(5, 35)
Model: y = a + b e^(x/4)
a: 3.056221
b: 1.922654
Curve Equation: y = 3.056221 + 1.922654 e^(x/4)

```

---

## Numerical Integration

---

### Simpson's 1/3 Rule

#### Simpson 1/3 Theory
 
Simpson’s 1/3 Rule is a numerical method used to approximate the definite integral of a function. It estimates the area under a curve by approximating the function with quadratic polynomials over subintervals. It is more accurate than the Trapezoidal Rule for the same number of intervals.
  
- Requires the number of subintervals `n` to be **even**.  
- Uses **parabolic arcs** to approximate the curve of the function.  
- Step size is calculated as `h = (b - a)/n`.  
- Formula combines function values at endpoints and midpoints of subintervals.

**Formula:**  
∫[a to b] f(x) dx ≈ (h/3) * [f(x0) + 4(f(x1) + f(x3) + ... + f(xn-1)) + 2(f(x2) + f(x4) + ... + f(xn-2)) + f(xn)]


Where:  
- `x0 = a`, `xn = b`  
- `h` is the step size  
- `n` is even

**Example:**  
Approximate ∫0^2 (1 + x^2) dx with n = 4:  
1. Step size: `h = (2 - 0)/4 = 0.5`  
2. Points: `x0=0, x1=0.5, x2=1, x3=1.5, x4=2`  
3. Apply formula:  
∫0^2 (1 + x^2) dx ≈ (0.5/3) * [f(0) + 4*(f(0.5) + f(1.5)) + 2*f(1) + f(2)]




#### Simpson 1/3 Code

```cpp
#include<bits/stdc++.h>
using namespace std;

typedef double db;

db poly(db x, vector<db>& c){
    db res=0, p=1;
    for(int i=0;i<c.size();i++){
        res+=c[i]*p;
        p*=x;
    }
    return res;
}

int main(){
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int t;
    fin>>t;
    fout<<"Total Test Cases: "<<t<<"\n\n";

    for(int tc=1;tc<=t;tc++){
        int deg;
        fin>>deg;

        vector<db> coef(deg+1);
        for(int i=0;i<=deg;i++) fin>>coef[i];

        db a,b;
        int n;
        fin>>a>>b>>n;

        fout<<"Test Case "<<tc<<"\n";

        if(n%2!=0){
            fout<<"Invalid input: n must be even\n\n";
            continue;
        }

        db h=(b-a)/n;
        db sum=poly(a,coef)+poly(b,coef);

        for(int i=1;i<n;i++){
            db x=a+i*h;
            if(i%2==0) sum+=2*poly(x,coef);
            else sum+=4*poly(x,coef);
        }

        db ans=(h/3)*sum;


        if(tc%3==1){
            fout<<"Integration Interval: ["<<a<<", "<<b<<"]\n";
            fout<<"Polynomial Order: "<<deg<<"\n";
            fout<<"Given Coefficients: ";
            for(db v:coef) fout<<v<<" ";
            fout<<"\nTotal Segments Used: "<<n<<"\n";
            fout<<"Computed Step Size h: "<<fixed<<setprecision(6)<<h<<"\n";
            fout<<"Final Integral Approximation: "<<fixed<<setprecision(6)<<ans<<"\n\n";
        }
        else if(tc%3==2){
            fout<<"Polynomial Order: "<<deg<<"\n";
            fout<<"Integration Interval: ["<<a<<", "<<b<<"]\n";
            fout<<"Total Segments Used: "<<n<<"\n";
            fout<<"Given Coefficients: ";
            for(db v:coef) fout<<v<<" ";
            fout<<"\nComputed Step Size h: "<<fixed<<setprecision(6)<<h<<"\n";
            fout<<"Final Integral Approximation: "<<fixed<<setprecision(6)<<ans<<"\n\n";
        }
        else{
            fout<<"Given Coefficients: ";
            for(db v:coef) fout<<v<<" ";
            fout<<"\nPolynomial Order: "<<deg<<"\n";
            fout<<"Integration Interval: ["<<a<<", "<<b<<"]\n";
            fout<<"Computed Step Size h: "<<fixed<<setprecision(6)<<h<<"\n";
            fout<<"Total Segments Used: "<<n<<"\n";
            fout<<"Final Integral Approximation: "<<fixed<<setprecision(6)<<ans<<"\n\n";
        }
    }

    fin.close();
    fout.close();
    return 0;
}


```

#### Simpson 1/3 Input

```
3
2
1 0 1
0 2 4
3
0 2 0 1
1 3 6
2
2 -1 1
0 2 8


```

#### Simpson 1/3 Output

```
Total Test Cases: 3

Test Case 1
Integration Interval: [0, 2]
Polynomial Order: 2
Given Coefficients: 1 0 1 
Total Segments Used: 4
Computed Step Size h: 0.500000
Final Integral Approximation: 4.666667

Test Case 2
Polynomial Order: 3
Integration Interval: [1.000000, 3.000000]
Total Segments Used: 6
Given Coefficients: 0.000000 2.000000 0.000000 1.000000 
Computed Step Size h: 0.333333
Final Integral Approximation: 28.000000

Test Case 3
Given Coefficients: 2.000000 -1.000000 1.000000 
Polynomial Order: 2
Integration Interval: [0.000000, 2.000000]
Computed Step Size h: 0.250000
Total Segments Used: 8
Final Integral Approximation: 4.666667



```

---

### Simpson's 3/8 Rule

#### Simpson 3/8 Theory
 
Simpson’s 3/8 Rule is a numerical method used to approximate the definite integral of a function. It is similar to Simpson’s 1/3 Rule but uses cubic polynomials to approximate the curve over subintervals. It is useful when the number of subintervals is a multiple of 3.
  
- Requires the number of subintervals `n` to be a **multiple of 3**.  
- Uses **cubic polynomials** to approximate the curve.  
- Step size is calculated as `h = (b - a)/n`.  
- More accurate than the Trapezoidal Rule for smooth functions.

**Formula:**  
∫[a to b] f(x) dx ≈ (3h/8) * [f(x0) + 3(f(x1) + f(x2) + f(x4) + f(x5) + ...) + 2(f(x3) + f(x6) + ...) + f(xn)]

Where:  
- `x0 = a`, `xn = b`  
- `h` is the step size  
- `n` is a multiple of 3

**Example:**  
Approximate ∫0^3 (1 + x^2) dx with n = 3:  
1. Step size: `h = (3 - 0)/3 = 1`  
2. Points: `x0=0, x1=1, x2=2, x3=3`  
3. Apply formula:  
∫0^3 (1 + x^2) dx ≈ (31/8) * [f(0) + 3(f(1) + f(2)) + f(3)]



#### Simpson 3/8 Code

```cpp
#include<bits/stdc++.h>
using namespace std;

typedef double db;


db poly(db x, vector<db>& coeff){
    db res = 0;
    db p = 1;
    for(int i=0;i<coeff.size();i++){
        res += coeff[i]*p;
        p *= x;
    }
    return res;
}

int main(){
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int t;
    fin >> t;
    fout << "Total Test Cases: " << t << "\n\n";

    for(int tc=1; tc<=t; tc++){
        int deg;
        fin >> deg;

        vector<db> coeff(deg+1);
        for(int i=0;i<=deg;i++) fin >> coeff[i];

        db a,b;
        int n;
        fin >> a >> b >> n;

        fout << "Test Case #" << tc << "\n";

        if(n%3 != 0){
            fout << "Invalid input: number of intervals must be multiple of 3\n\n";
            continue;
        }

        db h = (b-a)/n;
        db sum = poly(a, coeff) + poly(b, coeff);

        for(int i=1; i<n; i++){
            db x = a + i*h;
            if(i%3 == 0) sum += 2*poly(x, coeff);
            else sum += 3*poly(x, coeff);
        }

        db result = (3*h/8)*sum;


        fout << "Polynomial degree: " << deg << "\n";
        fout << "Coefficients: ";
        for(db c : coeff) fout << c << " ";
        fout << "\nInterval: [" << a << ", " << b << "]\n";
        fout << "Number of intervals: " << n << "\n";
        fout << "Step size (h): " << fixed << setprecision(6) << h << "\n";
        fout << "Integral result: " << fixed << setprecision(6) << result << "\n\n";
    }

    fin.close();
    fout.close();
    return 0;
}

```

#### Simpson 3/8 Input

```
3
2
1 0 1
0 2 6
3
0 1 0 0
1 3 9
2
2 -1 1
0 2 12

```

#### Simpson 3/8 Output

```
Total Test Cases: 3

Test Case #1
Polynomial degree: 2
Coefficients: 1 0 1 
Interval: [0, 2]
Number of intervals: 6
Step size (h): 0.333333
Integral result: 4.666667

Test Case #2
Polynomial degree: 3
Coefficients: 0.000000 1.000000 0.000000 0.000000 
Interval: [1.000000, 3.000000]
Number of intervals: 9
Step size (h): 0.222222
Integral result: 4.000000

Test Case #3
Polynomial degree: 2
Coefficients: 2.000000 -1.000000 1.000000 
Interval: [0.000000, 2.000000]
Number of intervals: 12
Step size (h): 0.166667
Integral result: 4.666667


```


