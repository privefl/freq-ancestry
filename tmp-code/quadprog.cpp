//  THIS IS BASED ON THE FORTRAN CODE OF R PACKAGE QUADPROG, REWRITTEN IN RCPP
//  IT ALSO ASSUMES THAT DMAT IS ALREADY FACTORIZED
//
//  Rcpp Port
//  Copyright (C) 2024 Florian Privé <florian.prive.21@gmail.com>
//  Original Fortran Code
//  Copyright (C) 1995-2010 Berwin A. Turlach <Berwin.Turlach@gmail.com>
//
//  this routine uses the Goldfarb/Idnani algorithm to solve the
//  following minimization problem:
//
//        minimize  -d^T x + 1/2 *  x^T D x
//        where   A1^T x  = b1
//                A2^T x >= b2
//
//  the matrix D is assumed to be positive definite.  Especially,
//  w.l.o.g. D is assumed to be symmetric.
//
//  Input parameter:
//  dmat   nxn matrix, the matrix D from above (dp)
//         actually pass R^-1, where D=R^TR.
//         *** WILL BE DESTROYED ON EXIT ***
//  dvec   nx1 vector, the vector d from above (dp)
//         *** WILL BE DESTROYED ON EXIT ***
//         contains on exit the solution to the initial, i.e.,
//         unconstrained problem
//  amat   nxq matrix, the matrix A from above (dp) [ A=(A1 A2)^T ]
//         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
//             CHANGED SIGNES ON EXIT ***
//  bvec   qx1 vector, the vector of constants b in the constraints (dp)
//         [ b = (b1^T b2^T)^T ]
//         *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
//             CHANGED SIGNES ON EXIT ***
//  meq    integer, the number of equality constraints, 0 <= meq <= q.
//
//
List qpgen2(NumericMatrix& dmat,
            NumericVector& dvec,
            NumericMatrix& amat,
            NumericVector& bvec,
            int meq) {  // TODO: -1 WHEN PASSING TO C++

  int n = dmat.nrow(), q = amat.ncol();
  int r = std::min(n, q);
  int l = 2 * n + r * (r + 5) / 2 + 2 * q + 1;

  // double gc, gs, nu

  IntegerVector iact(q);  // the constraints which are active in the final fit
  NumericVector sol(n);   // the final solution (x in the notation above)
  NumericVector lagr(q);  // the final Lagrange multipliers
  NumericVector work(l);  // working vector with length of at least l

  // code gleaned from Powell's ZQPCVX routine to determine a small number that can be
  // assumed to be an upper bound on the relative precision of the computer arithmetic.
  // (not sure it is still required though)
  double vsmall = 1.0e-60; while ((1.0 + 0.1 * vsmall) <= 1.0) vsmall += vsmall;

  // store the initial dvec to calculate below the unconstrained minima of the critical value.
  for (int i = 0; i < n; i++) work[i] = dvec[i];

  // get the initial solution (multiply d first with R^-T and then with R^-1)
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      sol[j] += dmat(i, j) * dvec[i];

  for (int j = 0; j < n; j++) {
    dvec[j] = 0;
    for (int i = 0; i < n; i++)
      dvec[j] += dmat(j, i) * sol[i];
  }

  // set lower triangular of dmat to zero, store dvec in sol and
  // calculate value of the criterion at unconstrained minima
  double crval = 0;       // the value of the criterion at the minimum
  for (int j = 0; j < n; j++) {
    sol[j] = dvec[j];
    crval += work[j] * sol[j];
    work[j] = 0;
    for (int i = j + 1; i < n; i++) dmat(i, j) = 0;
  }
  crval *= -0.5;

  // calculate some constants, i.e., from which index on the different
  // quantities are stored in the work matrix
  int iwzv  = n;
  int iwrv  = iwzv + n;
  int iwuv  = iwrv + r;
  int iwrm  = iwuv + r + 1;
  int iwsv  = iwrm + r * (r + 1) / 2;
  int iwnbv = iwsv + q;

  // calculate the norm of each column of the A matrix
  for (int i = 0; i < q; i++) {
    double sum = 0;
    for (int j = 0; j < n; j++) sum += amat(j, i) * amat(j, i);
    work[iwnbv + i] = ::sqrt(sum);
  }

  int nact = 0;           // the number of constraints active in the final fit
  int iter_1 = 0;         // the number of "main" iterations
  int iter_2 = 0;         // how many constraints were deleted after they became active

  while (true) {  // loop 50

    // start a new iteration
    iter_1++;

    // calculate all constraints and check which are still violated
    // for the equality constraints we have to check whether the normal
    // vector has to be negated (as well as bvec in that case)
    for (int i = 0; i < q; i++) {

      double sum = -bvec[i];
      for (int j = 0; j < n; j++) sum += amat(j, i) * sol[j];
      if (sum < vsmall) sum = 0;

      int l = iwsv + i;
      if (i > meq) {
        work[l] = sum;
      } else {
        work[l] = -std::abs(sum);
        if (sum > 0) {
          // switch sign of bvec[i] and amat.col(i)
          bvec[i] *= -1;
          for (int j = 0; j < n; j++) amat(j, i) *= -1;
        }
      }
    }

    // as safeguard against rounding errors, set already active constraints explicitly to 0
    for (int i = 0; i < nact; i++) work[iwsv + iact[i]] = 0;

    // we weight each violation by the number of non-zero elements in the
    // corresponding row of A. then we choose the violated constraint which
    // has maximal absolute value, i.e., the minimum.
    // by obvious commenting and uncommenting we can choose the strategy to
    // take always the first constraint which is violated. ;-)
    int nvl = -1;
    double temp = 0;
    for (int i = 0; i < q; i++) {
      if (work[iwsv + i] < temp * work[iwnbv + i]) {
        nvl = i;
        temp = work[iwsv + i] / work[iwnbv + i];
      }
      // if (work[iwsv + i] < 0) {
      //   nvl = i;
      //   break;
      // }
    }
    if (nvl == -1) {
      for (int i = 0; i < nact; i++) lagr[iact[i]] = work[iwuv + i];
      break;  // TODO: RETURN?
    }

    // calculate d=J^Tn^+ where n^+ is the normal vector of the violated
    // constraint. J is stored in dmat in this implementation!!
    // if we drop a constraint, we have to jump back here.
    while (true) {  // loop 55

      for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < n; j++) sum += dmat(j, i) * amat(j, nvl);
        work[i] = sum;
      }

      // now calculate z = J_2 d_2
      int l1 = iwzv;
      for (int i = 0; i < n; i++) work[l1 + i] = 0;
      for (int j = nact + 1; j < n; j++)
        for (int i = 0; i < n; i++)
          work[l1 + i] += dmat(i, j) * work[j];

      // and r = R^{-1} d_1, check also if r has positive elements
      // (among the entries corresponding to inequalities constraints).
      bool t1inf = true;
      double t1;          // the min of u/r
      int it1;            // index in which component t1 occurs
      for (int i = nact - 1; i >= 0; i--) {
        double sum = work[i];
        int l = iwrm + i * (i + 3) / 2;
        int l1 = l - i;
        for (int j = i + 1; j < nact; j++) {
          sum -= work[l] * work[iwrv + j];
          l += j;
        }
        sum /= work[l1];
        work[iwrv + i] = sum;
        if (iact[i] > meq && sum > 0) {
          t1inf = false;
          it1 = i;
        }
      }

      // if r has positive elements, find the partial step length t1, which is
      // the maximum step in dual space without violating dual feasibility.
      if (!t1inf) {
        t1 = work[iwuv + it1] / work[iwrv + it1];
        for (int i = 0; i < nact; i++) {
          if (iact[i] > meq && work[iwrv + i] > 0) {
            temp = work[iwuv + i] / work[iwrv + i];
            if (temp < t1) {
              t1 = temp;
              it1 = i;
            }
          }
        }
      }

      // test if the z vector is equal to zero
      double sum = 0;
      for (int i = 0; i < n; i++) {
        int i2 = iwzv + i;
        sum += work[i2] * work[i2];
      }
      if (sum <= vsmall) {
        // No step in primal space such that the new constraint becomes feasible.
        if (t1inf) Rcpp::stop("constraints are inconsistent, no solution!");

        // we take a partial step in dual space and drop constraint it1,
        // that is, we drop the it1-th active constraint.
        // then we continue at step 2(a) (marked by label 55)
        for (int i = 0; i < nact; i++) work[iwuv + i] -= t1 * work[iwrv + i];
        work[iwuv + nact + 1] += t1;
        break;  // of loop 55
      } else {
        // compute full step length t2, minimum step in primal space such that the
        // constraint becomes feasible. keep sum (which is z^Tn^+) to update crval below!
        double sum = 0;
        for (int i = 0; i < n; i++) sum += work[iwzv + i] * amat(i, nvl);
        double tt = -work[iwsv + nvl] / sum;
        bool t2min = true;
        if (!t1inf && t1 < tt) {
          tt = t1;
          t2min = false;
        }

        // take step in primal and dual space
        for (int i = 0; i < n; i++) sol[i] += tt * work[iwzv + i];
        crval += tt * sum * (tt / 2 + work[iwuv + nact + 1]);
        for (int i = 0; i < n; i++) work[iwuv + i] -= tt * work[iwrv + i];
        work[iwuv + nact + 1] += tt;

        // if it was a full step, then we check wheter further constraints are violated,
        // otherwise we can drop the current constraint and iterate once more
        if (t2min) {

          // we took a full step. Thus add constraint nvl to the list of active
          // constraints and update J and R
          nact++;
          iact[nact] = nvl;

          // to update R we have to put the first nact-1 components of the d vector
          // into column (nact) of R
          l = iwrm + (nact - 1) * nact / 2 + 1;
          for (int i = 0; i < (nact - 1); i++, l++) work[l] = work[i];

          // if now nact=n, then we just have to add the last element to the new
          // row of R.
          // Otherwise we use Givens transformations to turn the vector d(nact:n)
          // into a multiple of the first unit vector. That multiple goes into the
          // last element of the new row of R and J is accordingly updated by the
          // Givens transformations.
          if (nact == n) {
            work[l] = work[n];  // TODO: PROBLEM?
          } else {
            for (int i = n - 1; i >= nact, i--) {
              // we have to find the Givens rotation which will reduce the element
              // (l1) of d to zero.
              // if it is already zero we don't have to do anything, except of
              // decreasing l1
              if (work[i] == 0) continue;
              double gc = std::max(std::abs(work[i - 1]), std::abs(work[i]));
              double gs = std::min(std::abs(work[i - 1]), std::abs(work[i]));
              temp = std::abs(gc * ::sqrt(1 + gs / gc * gs gc));
              temp = (work[i - 1] >= 0) ? temp : -temp;
              gc = work[i - 1] / temp;
              gs = work[i]     / temp;

              // The Givens rotation is done with the matrix (gc gs, gs -gc).
              // If gc is one, then element (i) of d is zero compared with element
              // (l1-1). Hence we don't have to do anything.
              // If gc is zero, then we just have to switch column (i) and column (i-1)
              // of J. Since we only switch columns in J, we have to be careful how we
              // update d depending on the sign of gs.
              // Otherwise we have to apply the Givens rotation to these columns.
              // The i-1 element of d has to be updated to temp.
              if (gc == 1) continue;
              if (gc == 0) {
                work[i - 1] = gs * temp;
                for (int j = 0; j < n; j++) {
                  temp           = dmat(j, i - 1);
                  dmat(j, i - 1) = dmat(j, i);
                  dmat(j, i)     = temp;
                }
              } else {
                work[i - 1] = temp;
                double nu = gs / (1 + gc);
                for (int j = 0; j < n; j++) {
                  temp           = gc * dmat(j, i - 1) + gs * dmat(j, i);
                  dmat(j, i)     = nu * (dmat(j, i - 1) + temp) - dmat(j, i);
                  dmat(j, i - 1) = temp;
                }
              }
            }

            // l is still pointing to element (nact,nact) of the matrix R.
            // So store d(nact) in R(nact,nact)
            work[l] = work[nact];
          }

        } else {  // not t2min

          // we took a partial step in dual space. Thus drop constraint it1,
          // that is, we drop the it1-th active constraint.
          // then we continue at step 2(a) (marked by label 55)
          // but since the fit changed, we now have to recalculate "how much"
          // the fit violates the chosen constraint.
          double sum = -bvec[nvl];
          for (int j = 0; j < n; j++) sum += sol[j] * amat(j, nvl);
          if (nvl > meq) {
            work[iwsv + nvl] = sum;
          } else {
            work[iwsv + nvl] = -std::abs(sum);
            if (sum > 0) {
              // switch sign of bvec[nvl] and amat.col(nvl)
              bvec[nvl] *= -1;
              for (int j = 0; j < n; j++) amat(j, nvl) *= -1;
            }
          }
          break;  // of loop 55
        }
      }

    }
    //
    // Drop constraint it1
    //
    700  continue
    //
    // if it1 = nact it is only necessary to update the vector u and nact
    //
    if (it1 .EQ. nact) goto 799
    //
    // After updating one row of R (column of J) we will also come back here
    //
    797  continue
    //
    // we have to find the Givens rotation which will reduce the element
    // (it1+1,it1+1) of R to zero.
    // if it is already zero we don't have to do anything except of updating
    // u, iact, and shifting column (it1+1) of R to column (it1)
    // l  will point to element (1,it1+1) of R
    // l1 will point to element (it1+1,it1+1) of R
    //
    l  = iwrm + (it1*(it1+1))/2 + 1
    l1 = l+it1
    if (work(l1) .EQ. 0.d0) goto 798
    gc   = max(abs(work(l1-1)),abs(work(l1)))
      gs   = min(abs(work(l1-1)),abs(work(l1)))
      temp = sign(gc*sqrt(1+(gs/gc)*(gs/gc)), work(l1-1))
      gc   = work(l1-1)/temp
    gs   = work(l1)/temp
    //
    // The Givens rotatin is done with the matrix (gc gs, gs -gc).
    // If gc is one, then element (it1+1,it1+1) of R is zero compared with
    // element (it1,it1+1). Hence we don't have to do anything.
    // if gc is zero, then we just have to switch row (it1) and row (it1+1)
    // of R and column (it1) and column (it1+1) of J. Since we swithc rows in
    // R and columns in J, we can ignore the sign of gs.
    // Otherwise we have to apply the Givens rotation to these rows/columns.
    //
    if (gc .EQ. 1.d0) goto 798
    if (gc .EQ. 0.d0) then
    do 710 i=it1+1,nact
    temp       = work(l1-1)
      work(l1-1) = work(l1)
      work(l1)   = temp
    l1 = l1+i
    710     continue
    do 711 i=1,n
    temp          = dmat(i,it1)
      dmat(i,it1)   = dmat(i,it1+1)
      dmat(i,it1+1) = temp
    711     continue
    else
      nu = gs/(1.d0+gc)
      do 720 i=it1+1,nact
      temp       = gc*work(l1-1) + gs*work(l1)
        work(l1)   = nu*(work(l1-1)+temp) - work(l1)
        work(l1-1) = temp
      l1 = l1+i
      720     continue
      do 721 i=1,n
      temp          = gc*dmat(i,it1) + gs*dmat(i,it1+1)
        dmat(i,it1+1) = nu*(dmat(i,it1)+temp) - dmat(i,it1+1)
        dmat(i,it1)   = temp
      721     continue
      endif
      //
      // shift column (it1+1) of R to column (it1) (that is, the first it1
      // elements). The posit1on of element (1,it1+1) of R was calculated above
      // and stored in l.
      //
      798  continue
      l1 = l-it1
      do 730 i=1,it1
      work(l1)=work(l)
        l  = l+1
      l1 = l1+1
      730  continue
      //
      // update vector u and iact as necessary
      // Continue with updating the matrices J and R
      //
      work(iwuv+it1) = work(iwuv+it1+1)
        iact(it1)      = iact(it1+1)
        it1 = it1+1
      if (it1 .LT. nact) goto 797
      799  work(iwuv+nact)   = work(iwuv+nact+1)
        work(iwuv+nact+1) = 0.d0
      iact(nact)        = 0
      nact = nact-1
      iter_2++;
      goto 55
  }
  999  continue
  return
  end

  list(solution=res1$sol,
       value=res1$crval,
       unconstrained.solution=res1$dvec,
       iterations=res1$iter,
       Lagrangian = res1$lagr,
       iact=res1$iact[1:res1$nact])


}

/*** R
timesTwo(42)
*/
