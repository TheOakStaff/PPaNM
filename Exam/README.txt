--------------------------------------------------------------------------------
      As my student number ends in 08 I have made exam question number 8.
--------------------------------------------------------------------------------

  I have made an implementation of the one-sided Jacobi algorithm for single
value decomposition algorithm. The algorithm takes any given square-real matrix
'A' and factorizes it into two unitarian matrices U and V (U^TU = UU^T = I,
V^TV = VV^T = I), as well as a diagonal matrix wit the singular values of
matrix A.
  The resulting matrices (A',U,V,D) of my algorithm as well as the checks for
A = UDV^T, U^TU = UU^T = I & V^TV = VV^T = I can be found in the out.txt file.

  Additionally, I have compared the speed of my implementation of the one-sided
Jacobi algorithm for single value decomposition with the one provided by the
GSL - GNU Scientific Library (v. 2.5) "gsl_linalg_SV_decomp_jacobi". Both
algorithms are given the same initial matrix of size n and the time it takes to
complete the operation is then measured. Generally results indicate that the two
algorithms are evenly matched for small matrices n < 5, for larger matrices the
GSL algorithm is noticeably faster, at n = 400 the difference is one order of
magnitude.

I rate my algorithm a 8/10.

--------------------------------------------------------------------------------
                              Exam question 8
--------------------------------------------------------------------------------

-- Title --
One-sided Jacobi algorithm for Singular Value Decomposition

-- Introduction --
The singular value decomposition (SVD) of a (real square, for simplicity) matrix
A is a representation of the matrix in the form

A = U D V^T

where matrix D is diagonal with non-negative elements and matrices U and V are
orghogonal. The diagonal elements of matrix D can always be chosen non-negative
by multiplying the relevant columns of matrix U with (-1).

SVD can be used to solve a number of problems in linear algebra.

-- Problem --
Implement the one-sided Jacobi SVD algorithm.

-- Algorithm --
In this method the elementary iteration is given as
A → A J(θ,p,q)

where the indices (p,q) are swept cyclicly (p=1..n, q=p+1..n) and where the
angle θ is chosen such that the columns number p and q of the matrix AJ(θ,p,q)
are orthogonal. One can show that the angle should be taken from the following
equation (you should use atan2 function),
tan(2θ)=2ap^Taq /(aq^Taq - ap^Tap)

where ai is the i-th column of matrix A (check this).
After the iterations converge and the matrix A'=AJ (where J is the accumulation
of the individual rotations) has orthogonal columns, the SVD is simply given as

A=U D V^T

where
V=J, Dii=||a'i||, ui=a'i/||a'i||,

where a'i is the i-th column of matrix A' and ui us the i-th column of matrix U.
--------------------------------------------------------------------------------
