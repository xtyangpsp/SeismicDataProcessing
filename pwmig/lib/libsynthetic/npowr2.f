      function npowr2(n)
c
c finds the next power of 2 .ge.n
c
      ipowr=alog10(2.*float(n)-1.)/.301029996
      if(n.eq.1) ipowr=1
c--changed by glp feb 2013 as this requires a fortran library
c--routine that complicates linking to C
c     npowr2=2**ipowr
      npowr2=1;
      do 100 i=1,ipowr;
           npowr2=2*npowr2
  100 continue
      return
      end
