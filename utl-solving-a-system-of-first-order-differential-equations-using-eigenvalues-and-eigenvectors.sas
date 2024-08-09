%let pgm=utl-solving-a-system-of-first-order-differential-equations-using-eigenvalues-and-eigenvectors;

Solving a system of first order differential equations using eigenvalues and eigenvectors

Please be skeptical of all the solutions I propose in github, github is a development site
and welcomes issues and corrections. I make mistakes.

   1 symbolic solution sympy
   2 r numerical solution
   3 comparison sympy & r
   4 sas normal non-normal eigen data

github
https://github.com/rogerjdeangelis/utl-solving-a-system-of-first-order-differential-equations-using-eigenvalues-and-eigenvectors
https://github.com/rogerjdeangelis/utl-solving-a-system-of-first-order-differential-equations-using-eigenvalues-and-eigenvectors

/*     _     _
| |__ (_)___| |_ ___  _ __ _   _
| `_ \| / __| __/ _ \| `__| | | |
| | | | \__ \ || (_) | |  | |_| |
|_| |_|_|___/\__\___/|_|   \__, |
                           |___/
*/

Early Developments and Influences of Symbolic Mathematics

Before Mathematica, there were several research systems aimed at performing algebraic
computations using computers. Stephen Wolfram, who started his career as a physicist,
found these systems inadequate for his needs in theoretical physics.
This led him to develop his own system in the late 1970s,
which eventually became the first commercial computer algebra system.

The development of symbolic mathematics also has philosophical roots,
as explored in the works of Edmund Husserl and Jacob Klein.
These philosophers examined the nature of symbolic numbers and their ontological
assumptions, which influenced the modern understanding of symbolic mathematics.
Klein, for instance, analyzed the evolution of symbolic
algebra from its origins in Diophantine algebra to its modern form,
which has implications for mathematical physics.

Diophantine algebra refers to the study of Diophantine equations,
which are polynomial equations that involve two or more unknowns
with integer coefficients

Related repo (Diophantine equations?)
https://tinyurl.com/ytnv4jk7
https://github.com/rogerjdeangelis/utl-fun-with-sympy-infinite-series-and-integrals-to-define-common-functions-and-constants

/*               _     _
 _ __  _ __ ___ | |__ | | ___ _ __ ___
| `_ \| `__/ _ \| `_ \| |/ _ \ `_ ` _ \
| |_) | | | (_) | |_) | |  __/ | | | | |
| .__/|_|  \___/|_.__/|_|\___|_| |_| |_|
|_|
*/

/**************************************************************************************************************************/
/*                                                                                                                        */
/*   Solve for x(t) and y(t)                                                                                              */
/*                                                                                                                        */
/*   Consider a system of linear differential equations:                                                                  */
/*                                                                                                                        */
/*          dx                                                                                                            */
/*          --     =  a*x  + b*y                                                                                          */
/*          dt                                                                                                            */
/*                                                                                                                        */
/*          dy                                                                                                            */
/*          --     =  b*x  + a*y                                                                                          */
/*          dt                                                                                                            */
/*                                                                                                                        */
/*   Suppose                                                                                                              */
/*                                                                                                                        */
/*      dx/dt measures the x component of the velocity of a billard ball                                                  */
/*      dy/dt measures the y component of the velocity of a billard ball                                                  */
/*                                                                                                                        */
/*      y(t) is the y position                                                                                            */
/*      x(t) is the x position                                                                                            */
/*                                                                                                                        */
/*   For the R we to assume a solution involving  exp(rx) because                                                         */
/*   the intergral and derivative of exp(rx) are functions of exp(rx).                                                    */
/*   This assumptions helps solve a wide variatey of differential equations                                               */
/*                                                                                                                        */
/*   We do nit need to make the assumption when using sympy                                                               */
/*   Sympy has built in knowlege of techniques to solve dirrential equations                                              */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*                       _           _ _                  _       _   _
/ |  ___ _   _ _ __ ___ | |__   ___ | (_) ___   ___  ___ | |_   _| |_(_) ___  _ __
| | / __| | | | `_ ` _ \| `_ \ / _ \| | |/ __| / __|/ _ \| | | | | __| |/ _ \| `_ \
| | \__ \ |_| | | | | | | |_) | (_) | | | (__  \__ \ (_) | | |_| | |_| | (_) | | | |
|_| |___/\__, |_| |_| |_|_.__/ \___/|_|_|\___| |___/\___/|_|\__,_|\__|_|\___/|_| |_|
         |___/
*/


%utl_pybeginx;
parmcards4;
import sympy as sp
from sympy import symbols, Function, Eq, dsolve, Derivative, pprint, dsolve

# Define the variables and functions
t = sp.symbols('t')
x = sp.Function('x')(t)
y = sp.Function('y')(t)

# Define the coefficients of the system
a, b, c, d = sp.symbols('a b c d')

# Define the system matrix
A = sp.Matrix([[a, b], [b, a]])
pprint(A)
# Find the eigenvalues and eigenvectors
eigenvals = A.eigenvals()
eigenvects = A.eigenvects()
pprint(eigenvects)

# Display eigenvalues and eigenvectors
print("Eigenvalues:")
for val, mult in eigenvals.items():
    pprint(f"Eigenvalue: {val}, Multiplicity: {mult}")

print("\nEigenvectors:")
for val, mult, vects in eigenvects:
    pprint(f"Eigenvalue: {val}, Multiplicity: {mult}, Eigenvectors: {vects}")

# Solve the system using the eigenvectors and eigenvalues
C1, C2 = sp.symbols('C1 C2')
solutions = []

for val, mult, vects in eigenvects:
    for vect in vects:
        sol = sp.exp(val * t) * vect
        solutions.append(sol)

# General solution
general_solution = C1 * solutions[0] + C2 * solutions[1]
x_sol, y_sol = general_solution

print("\nGeneral solution:")
print(f"x(t) = {x_sol}")
print(f"y(t) = {y_sol}")
;;;;
%utl_pyendx;


/*           _               _
  ___  _   _| |_ _ __  _   _| |_
 / _ \| | | | __| `_ \| | | | __|
| (_) | |_| | |_| |_) | |_| | |_
 \___/ \__,_|\__| .__/ \__,_|\__|
                |_|
*/


/**************************************************************************************************************************/
/*                                                                                                                        */
/*  Eigen Values                                                                                                          */
/*                                                                                                                        */
/*        a-b  = 3                                                                                                        */
/*        a+b  = 7                                                                                                        */
/*                                                                                                                        */
/*  Raw non-normal eiegenvectors (because of linear dependence we use v11=1 and v12=1 to and solve for v12 and v22)       */
/*                                                                                                                        */
/*    --  --                                                                                                              */
/*   | 1 -1 1                                                                                                             */
/*   | 1  1 1                                                                                                             */
/*    --   --                                                                                                             */
/*                                                                                                                        */
/*               t*(a - b)        t*(a + b)                                                                               */
/*  x(t) = - C1*e           + C2*e                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/*               t*(a - b)        t*(a + b)                                                                               */
/*  y(t) =   C1*e           + C2*e                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/*  Initial conditions                                                                                                    */
/*                                                                                                                        */
/*  x(t=0) =  c1 + c2                                                                                                     */
/*  y(t=0) =  c1 + c2                                                                                                     */
/*                                                                                                                        */
/*  Here is one solution                                                                                                  */
/*                                                                                                                        */
/*  x(t) = -C1*exp(t*(a - b)) + C2*exp(t*(a + b))                                                                         */
/*                                                                                                                        */
/*  y(t) =  C1*exp(t*(a - b)) + C2*exp(t*(a + b))                                                                         */
/*                                                                                                                        */
/*                                                                                                                        */
/*  Note C1 and C2 are arbirary and are the constants of integration                                                      */
/*  and are arbitraty                                                                                                     */
/*                                                                                                                        */
/*                                                                                                                        */
/*  Lets  Evaluate at couple of points (c1 and c2 are arbitrary)                                                          */
/*                                                                                                                        */
/*    c1=1                                                                                                                */
/*    c2=1                                                                                                                */
/*                                                                                                                        */
/*     a=5                                                                                                                */
/*     b=2                                                                                                                */
/*     c=5                                                                                                                */
/*                                                                                                                        */
/*  x(0)   = -C1*exp(t*(a - b)) + C2*exp(t*(a + b)) = -c1 + c2 = 0                                                        */
/*  y(0)   =  C1*exp(t*(a - b)) + C2*exp(t*(a + b)) =  c1 = c2 = 2                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/*  x(-1)  = -C1*exp(t*(a - b)) + C2*exp(t*(a + b)) = -c1 + c2 =  -0.048875186                                            */
/*  y(-2)  =  C1*exp(t*(a - b)) + C2*exp(t*(a + b)) =  c1 = c2 =   0.05069895                                             */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*___                                           _           _             _       _   _
|___ \   _ __   _ __  _   _ _ __ ___   ___ _ __(_) ___ __ _| |  ___  ___ | |_   _| |_(_) ___  _ __
  __) | | `__| | `_ \| | | | `_ ` _ \ / _ \ `__| |/ __/ _` | | / __|/ _ \| | | | | __| |/ _ \| `_ \
 / __/  | |    | | | | |_| | | | | | |  __/ |  | | (_| (_| | | \__ \ (_) | | |_| | |_| | (_) | | | |
|_____| |_|    |_| |_|\__,_|_| |_| |_|\___|_|  |_|\___\__,_|_| |___/\___/|_|\__,_|\__|_|\___/|_| |_|

*/

%utl_rbeginx;
parmcards4;
# Define the matrix A
A <- matrix(c(5, 2, 2, 5), nrow = 2)

# Compute eigenvalues and eigenvectors
eigen_result <- eigen(A)
eigenvalues <- eigen_result$values
eigenvectors <- eigen_result$vectors
non_normalized_vectors <- apply(eigen_result$vectors, 2, function(v) v / v[2])
print(non_normalized_vectors)
# Display the results
cat("Eigenvalues:\n")
print(eigenvalues)
cat("non_normalized_vectors:\n")
print(non_normalized_vectors)

# Define the general solution
general_solution <- function(t, c1, c2) {
  v1 <- non_normalized_vectors[, 1]
  v2 <- non_normalized_vectors[, 2]
  lambda1 <- eigenvalues[1]
  lambda2 <- eigenvalues[2]
  y <- c1 * exp(lambda1 * t) * v1 + c2 * exp(lambda2 * t) * v2
  return(y)
}

# Example usage
c1 <- 1
c2 <- 1
solution <- general_solution(0, c1, c2)
cat("Solution at t = 0:\n")
print(solution)
solution <- general_solution(-1, c1, c2)
cat("Solution at t = -1:\n")
print(solution)
;;;;
%utl_rendx;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*  > print(non_normalized_vectors)                                                                                       */
/*       [,1] [,2]                                                                                                        */
/*  [1,]    1   -1                                                                                                        */
/*  [2,]    1    1                                                                                                        */
/*                                                                                                                        */
/*  > print(eigenvalues)                                                                                                  */
/*  [1] 7 3                                                                                                               */
/*                                                                                                                        */
/*  Solution at t = 0:                                                                                                    */
/*  > print(solution)                                                                                                     */
/*  [1] 0 2                                                                                                               */
/*  >                                                                                                                     */
/*  > print(solution)                                                                                                     */
/*  [1] -0.04887519  0.05069895                                                                                           */
/*                                                                                                                        */
/**************************************************************************************************************************/


/*____                                        _                                                     ___
|___ /    ___ ___  _ __ ___  _ __   __ _ _ __(_)___  ___  _ __   ___ _   _ _ __ ___  _ __  _   _   ( _ )    _ __
  |_ \   / __/ _ \| `_ ` _ \| `_ \ / _` | `__| / __|/ _ \| `_ \ / __| | | | `_ ` _ \| `_ \| | | |  / _ \/\ | `__|
 ___) | | (_| (_) | | | | | | |_) | (_| | |  | \__ \ (_) | | | |\__ \ |_| | | | | | | |_) | |_| | | (_>  < | |
|____/   \___\___/|_| |_| |_| .__/ \__,_|_|  |_|___/\___/|_| |_||___/\__, |_| |_| |_| .__/ \__, |  \___/\/ |_|
                            |_|                                      |___/          |_|    |___/

*/

Data check;

   a=5;   b=2;
   c=2;   d=5;

   t=-1;

   C1=1;
   C2=1;

   xt1  = -C1*exp(t*(a - b)) + C2*exp(t*(a + b));
   yt1  =  C1*exp(t*(a - b)) + C2*exp(t*(a + b));

   t=0;

   xt0  = -C1*exp(t*(a - b)) + C2*exp(t*(a + b));
   yt0  =  C1*exp(t*(a - b)) + C2*exp(t*(a + b));


   put 'sympy x(t=0) = ' xt0 / ' r numerical x(t=0) = 0' /;
   put 'sympy y(t=0) = ' yt0 / ' r numerical y(t=0) = 2' ;


   put  // 'sympy x(t=-1) = ' xt1 / ' r numerical x(t=0) = -0.04887519' /;
   put 'sympy y(t=-1) = ' yt1 / ' r numerical y(t=0) = 0.05069895' ;

run;quit;


/**************************************************************************************************************************/
/*                                                                                                                        */
/* R AGREES WITH SYMPY (Time=0)                                                                                           */
/* ============================                                                                                           */
/*                                                                                                                        */
/*  Time=0                                                                                                                */
/*                                                                                                                        */
/*   sympy x(t=0)       = 0                                                                                               */
/*   r numerical x(t=0) = 0                                                                                               */
/*                                                                                                                        */
/*   sympy y(t=0)       = 2                                                                                               */
/*   r numerical y(t=0) = 2                                                                                               */
/*                                                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/* R AGREES WITH SYMPY (Time=-1 yesterday)                                                                                */
/* ========================================                                                                               */
/*                                                                                                                        */
/*   sympy x(t=-1)      = -0.048875186                                                                                    */
/*   r numerical x(t=0) = -0.04887519                                                                                     */
/*                                                                                                                        */
/*   sympy y(t=-1)      = 0.0506989503                                                                                    */
/*   r numerical y(t=0) = 0.05069895                                                                                      */
/*                                                                                                                        */
/**************************************************************************************************************************/


/*  _                                                      _                                                            _
| || |    ___  __ _ ___   _ __   ___  _ __ _ __ ___   __ _| |  _ __   ___  _ __        _ __   ___  _ __ _ __ ___   __ _| |
| || |_  / __|/ _` / __| | `_ \ / _ \| `__| `_ ` _ \ / _` | | | `_ \ / _ \| `_ \ _____| `_ \ / _ \| `__| `_ ` _ \ / _` | |
|__   _| \__ \ (_| \__ \ | | | | (_) | |  | | | | | | (_| | | | | | | (_) | | | |_____| | | | (_) | |  | | | | | | (_| | |
   |_|   |___/\__,_|___/ |_| |_|\___/|_|  |_| |_| |_|\__,_|_| |_| |_|\___/|_| |_|     |_| |_|\___/|_|  |_| |_| |_|\__,_|_|
      _                        _       _
  ___(_) __ _  ___ _ __     __| | __ _| |_ __ _
 / _ \ |/ _` |/ _ \ `_ \   / _` |/ _` | __/ _` |
|  __/ | (_| |  __/ | | | | (_| | (_| | || (_| |
 \___|_|\__, |\___|_| |_|  \__,_|\__,_|\__\__,_|
        |___/
*/


 data _null_;;

    a=5;   b=2;
    c=2;   d=5;

    * eigenvalues;
    e1 = a/2 + d/2 - sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2;
    e2 = a/2 + d/2 + sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2;

    * non normalized eigenvectors;
    v11 =1;
    v21 =          -d/b + (a/2 + d/2 - sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2)/b;

    v12=1;
    v22 =          -d/b + (a/2 + d/2 + sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2)/b;

    * normalized eigervectors;
    norm11=v11   / sqrt(sum(1+v21**2));
    norm21=v21   / sqrt(sum(1+v21**2));
    norm12=v12   / sqrt(sum(1+v22**2));
    norm22=v22   / sqrt(sum(1+v22**2));

 put "Eigenvalues" / e1= / e2= //
     "Non-Normal eogennvectors" / v11= v21= / v12= v22= //
     "Normalized Eigenvectors" / norm11=  norm21= / norm12=  norm22= ;
 run;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/* Eigenvalues                                                                                                            */
/* e1=3                                                                                                                   */
/* e2=7                                                                                                                   */
/*                                                                                                                        */
/*                                                                                                                        */
/* I USED THESE IN SYMPY AND R                                                                                            */
/* ===========================                                                                                            */
/*                                                                                                                        */
/* Non-Normal eogennvectors                                                                                               */
/* v11=1 v21=-1                                                                                                           */
/* v12=1 v22=1                                                                                                            */
/*                                                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/* Normalized Eigenvectors                                                                                                */
/* norm11=0.7071067812 norm21=-0.707106781                                                                                */
/* norm12=0.7071067812 norm22=0.7071067812                                                                                */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*              _
  ___ _ __   __| |
 / _ \ `_ \ / _` |
|  __/ | | | (_| |
 \___|_| |_|\__,_|

*/
