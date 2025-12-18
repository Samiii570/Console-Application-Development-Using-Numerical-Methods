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

[Add your theory content here]

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

[Add your theory content here]

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

[Add your theory content here]

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

[Add your theory content here]

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
[Add your theory content here]

#### Bisection Code
```python
# Add your code here
````

#### Bisection Input

```
[Add your input format here]
```

#### Bisection Output

```
[Add your output format here]
```

---

### False Position Method

#### False Position Theory

[Add your theory content here]

#### False Position Code

```python
# Add your code here
```

#### False Position Input

```
[Add your input format here]
```

#### False Position Output

```
[Add your output format here]
```

---

### Newton Raphson Method

#### Newton Raphson Theory

[Add your theory content here]

#### Newton Raphson Code

```python
# Add your code here
```

#### Newton Raphson Input

```
[Add your input format here]
```

#### Newton Raphson Output

```
[Add your output format here]
```

---

### Secant Method

#### Secant Theory

[Add your theory content here]

#### Secant Code

```python
# Add your code here
```

#### Secant Input

```
[Add your input format here]
```

#### Secant Output

```
[Add your output format here]
```

---

## Differential Equation Solving

---

### Runge-Kutta 4th Order Method

#### Runge-Kutta Theory

[Add your theory content here]

#### Runge-Kutta Code

```python
# Add your code here
```

#### Runge-Kutta Input

```
[Add your input format here]
```

#### Runge-Kutta Output

```
[Add your output format here]
```

---

## Interpolation Methods

---

### Newton Forward Interpolation

#### Newton Forward Theory

[Add your theory content here]

#### Newton Forward Code

```python
# Add your code here
```

#### Newton Forward Input

```
[Add your input format here]
```

#### Newton Forward Output

```
[Add your output format here]
```

---

### Newton Backward Interpolation

#### Newton Backward Theory

[Add your theory content here]

#### Newton Backward Code

```python
# Add your code here
```

#### Newton Backward Input

```
[Add your input format here]
```

#### Newton Backward Output

```
[Add your output format here]
```

---

### Newton Divided Difference Interpolation

#### Newton Divided Difference Theory

[Add your theory content here]

#### Newton Divided Difference Code

```cpp

```

#### Newton Divided Difference Input

```
[Add your input format here]
```

#### Newton Divided Difference Output

```
[Add your output format here]
```

---

## Numerical Differentiation

---

### Differentiation by Forward Interpolation

#### Differentiation Forward Theory

[Add your theory content here]

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

[Add your theory content here]

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

[Add your theory content here]

#### Linear Regression Code

```python
# Add your code here
```

#### Linear Regression Input

```
[Add your input format here]
```

#### Linear Regression Output

```
[Add your output format here]
```

---

### Polynomial Regression

#### Polynomial Regression Theory

[Add your theory content here]

#### Polynomial Regression Code

```python
# Add your code here
```

#### Polynomial Regression Input

```
[Add your input format here]
```

#### Polynomial Regression Output

```
[Add your output format here]
```

---

### Transcendental Regression

#### Transcendental Regression Theory

[Add your theory content here]

#### Transcendental Regression Code

```python
# Add your code here
```

#### Transcendental Regression Input

```
[Add your input format here]
```

#### Transcendental Regression Output

```
[Add your output format here]
```

---

## Numerical Integration

---

### Simpson's 1/3 Rule

#### Simpson 1/3 Theory

```
Numerical integration is used to evaluate definite integrals when the analytical solution of an integral is either difficult or impossible to obtain. Simpson’s 1/3 Rule is one of the most widely used numerical integration techniques due to its simplicity and relatively high accuracy. In this method, the entire integration interval is divided into an even number of equal sub-intervals. The fundamental assumption of Simpson’s 1/3 Rule is that the function can be approximated by a second-degree polynomial (parabola) over each pair of consecutive sub-intervals.
Unlike the trapezoidal rule, which assumes linear variation of the function, Simpson’s 1/3 Rule assumes a smooth and continuous curvature of the function. As a result, it provides better accuracy, especially for functions that are reasonably smooth within the given limits.

Formula :
Let the lower and upper limits of integration be a and b respectively. The interval [a, b] is divided into n equal sub-intervals, where n must be an even number. 
step size, h = (b − a) / n.
The approximate value of the definite integral is given by:
∫a to b f(x) dx ≈ (h/3)[f(x0) + f(xn) + 4(f(x1)+f(x3)+…) + 2(f(x2)+f(x4)+…)]

Algorithm :
1. Read the lower limit a, upper limit b, and an even number of sub-intervals n. 
2. Calculate the step size h using h = (b − a) / n. 
3. Evaluate the function at the first and last points and add them.
4.Evaluate the function at all odd-indexed points and multiply each by 4. 
5. Evaluate the function at all even-indexed points (excluding boundaries) and multiply each by 2. 6. Add all the weighted function values. 
7. Multiply the total sum by h/3 to obtain the approximate value of the integral.

Working Principle :
Simpson’s 1/3 Rule works by fitting a parabolic curve through three equally spaced points of the function. By integrating this parabolic approximation instead of the original function, the area under the curve is estimated with improved precision. The use of different weighting factors for odd and even indexed points ensures that the curvature of the function is properly accounted for.

```

#### Simpson 1/3 Code

```cpp
#include<iostream>
#include<fstream>
#include<cmath>
using namespace std;

double f(double x){
    return x*x+1;
}

double simpsonOneThird(double a,double b,int n){
    double h=(b-a)/n;
    double sum=f(a)+f(b);
    for(int i=1;i<n;i++){
        double x=a+i*h;
        if(i%2==0) sum+=2*f(x);
        else sum+=4*f(x);
    }
    return (h/3)*sum;
}

int main(){
    ifstream fin("input.txt");
    ofstream fout("output_1_3.txt");
    double a,b;
    int n;
    int t=1;
    while(fin>>a>>b>>n){
        fout<<"Test case "<<t<<": a="<<a<<", b="<<b<<", n="<<n<<"\n";
        double result=simpsonOneThird(a,b,n);
        fout<<"Approximate integral="<<result<<"\n\n";
        t++;
    }
    fin.close();
    fout.close();
    return 0;
}


```

#### Simpson 1/3 Input

```
0 2 4
1 3 6
0 5 12

```

#### Simpson 1/3 Output

```
Test case 1: a=0, b=2, n=4
Approximate integral=4.66667

Test case 2: a=1, b=3, n=6
Approximate integral=10

Test case 3: a=0, b=5, n=12
Approximate integral=46.875

```

---

### Simpson's 3/8 Rule

#### Simpson 3/8 Theory

[Add your theory content here]

#### Simpson 3/8 Code

```cpp
# Add your code here
```

#### Simpson 3/8 Input

```
[Add your input format here]
```

#### Simpson 3/8 Output

```
[Add your output format here]
```


