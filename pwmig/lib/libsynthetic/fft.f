      subroutine fft(data,nn,isign)
c                                              a.shakal, 1/78, 10 jul 80
c        cooley-tukey 'fast fourier trnasform' in ansi fortran 77.
c
c           transform(j) = sum {data(i)*w**u(i-1)*(j-1)e}, where i and
c        j run from 1 to nn, and w = exp(sign*twopi*sqrtu-1e/nn).
c        data is a one-dimensional complex array (i.e., the real and
c        imaginary parts of the data are located immediately adjacent
c        in storage, such as fortran places them) whose length nn is
c        a power of two.  isign is +1 or -1, giving the sign of the
c        transform.  transform values are returned in array data,
c        replacing the input data.  the time is proportional to
c        n*log2(n), rather than the non-fft n**2.  modified from the
c        fortran ii coding from n.brenner's mit-ll tech rept.
c
      real data(1)
      pi = 3.1415926
c
      n = 2*nn
      j = 1
      do 5 i = 1,n,2
      if (.not.(i .lt. j)) go to 20001
      tempr = data(j)
      tempi = data(j+1)
      data(j) = data(i)
      data(j+1) = data(i+1)
      data(i) = tempr
      data(i+1) = tempi
20001 m = n/2
    3 if (.not.(j .gt. m)) go to 20004
      j = j-m
      m = m/2
      if (m .ge. 2) go to 3
20004 j = j+m
   5  continue
c
c
      mmax = 2
    6 if (.not.(mmax .ge. n)) go to 20007
      return
20007 if (.not.(mmax .lt. n)) go to 10001
      istep = 2*mmax
      pibymx = pi*float(isign)/float(mmax)
c
      do 8 m = 1,mmax,2
      theta = pibymx*float(m-1)
      wr = cos(theta)
      wi = sin(theta)
      do 8 i = m,n,istep
      j = i + mmax
      tempr = wr*data(j) - wi*data(j+1)
      tempi = wr*data(j+1) + wi*data(j)
      data(j) = data(i) - tempr
      data(j+1) = data(i+1) - tempi
      data(i) = data(i) + tempr
      data(i+1) = data(i+1) + tempi
   8  continue
      mmax = istep
      go to 6
10001 continue
20008 return
      end
