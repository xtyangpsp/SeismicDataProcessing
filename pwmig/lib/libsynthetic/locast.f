      subroutine locast(dirctn,invrse,forwrd)
      character dirctn*7
      logical forwrd,invrse
      integer ounit
      common /innout/ inunit,ounit
      if(dirctn.eq.'forward') go to 1
      if(dirctn.eq.'inverse') go to 2
      write(ounit,100)dirctn
  100 format(1x,a7,2x,'is meaningless to dfftr, use forward or inverse
     *only')
      invrse=.false.
      forwrd=.false.
      return
    1 invrse=.false.
      forwrd=.true.
      return
    2 invrse=.true.
      forwrd=.false.
      return
      end
