       subroutine kntsyn(nlyrs,alfm,betm,rhom,thikm,pr,fullrsp,
     *        dt,numpts,tsigma,wlevel,u0,w0,u1,w1,tn,
     *        rfr,ierr)
c--------------------------------------------------------------
c   interface routine derived from TJO program named respknt.  
c  This driver takes a layered earth model as input and returns
c  a set of vectors containing the synthetic seismogram for a plane wave
c  source (ray parameter p) illuminated from below (proto receiver
c  functions).  Note the vectors are fortran complex and if called
c  from C/C++ the caller needs to strip the real part stored in
c  alternate pairs (real is first of each pair).  
c
c  Original program was aimed at receiver functions and wrote output
c  sac files for vertical and radial.  This procedure returns the full
c  response set of 5 components defined by Kennett's original code.
c
c  arguments
c    nlyrs - number of layers (size of the next 4 vectors)
c    alfm,betm,rhom - P velocity, S velocity, and density respectively
c       (all of length nlyrs)
c    thikm - layer thickness in units of km.
c    pr - slowness of incident plane wave 
c    fullrsp - if nonzero compute all reverbs.  If zero primaries only
c    dt - sample interval (seconds)
c    numpts - number of points in output seismograms = size of arrays
c      that follow.  MUST be a power of 2 or the procedure will exit 
c      with ierr set nonzero (-2 actually)
c    tsigma - sigma of gaussian pulse in time domain.  Output is filtered
c      by gaussian with this width measure in time.
c    wlevel` - waterlevel deconvolution is used to produce the 
c      synthetic trace rfr using the vertical to deconvolve P to S radial.
c      This simulates standard data.  It is done here for efficiency
c      since synthetics are computed in the frequency domain..
c    u0,w0 - output synthetic seismograms for Z and R respectively incident P
c    u1,w1 - output synthetic seismogram for Z and R respectively incident S
c    tn - output synthetic for incident SH (transverse component)
c    ierr - error code (0 success, nonzero means something wrong)
c       -1 input model has more layers than internal work spaces allow
c       -2 numpts is not a power of 2 (required)
c
c  caution for naive user:  u0,w0,u1,w1, and tn must be allocated
c  externally.  Further, this procedure has the complex meaning 
c  the original array is 2*numpts in length.  When called from C/C++
c  caller must extract the real part from each vector as the first 
c  of each pair of the complex array.
c----------------------------------------------------------------
      implicit none
      real*4 tsigma
      real*4 t,dt,wlevel
      integer*4 nlyrs,numpts
      real*4 alfm(nlyrs),betm(nlyrs),rhom(nlyrs),thikm(nlyrs)
      real*4 pr
      integer*4 fullrsp
      complex u0(numpts),w0(numpts),u1(numpts),w1(numpts),tn(numpts)
      complex rfr(numpts)
      real*4 zpeak,zamp,rwlev
      real*4 fsigma
c--This is a time shift applied to rf estimate after decon.  Necessary
c--because spectral division cancels original phase.
      real*4 rftshift
      data rftshift/5.0/
      integer*4 maxlyr
      parameter(maxlyr=50)
      real qpm(maxlyr),qsm(maxlyr),ta(maxlyr),tb(maxlyr)
      integer*4 ierr
      complex dvp,dvs,drp,drs,dts,p,fr
      real freal
      real*8 wq,t1,t2,qa,qb,qabm,vabm
      real phdr
c  This monstrosity passes data to respknt procedure
      include 'kennet.inc'
      integer*2 cnv,rvb
      real*4 delf,fny
      integer*4 i,j,npts,nfpts,nft
c--new syntax needed for implicit none with this function
      integer :: npowr2
      real twopi
      data twopi/6.2831853/
      fsigma=1.0/tsigma;
      ierr=0
c--trap this error (added in conversion by glp from original)
      if(nlyrs.gt.maxlyr) then
            ierr=-1
            return
       endif

c--this appears to force stock values for attenuation parameters       
      do 1 i=1,nlyrs
       qpm(i) = 500.
       qsm(i) = 225.
       ta(i) = .16
       tb(i) = .26
 1    continue      
c
c     set up the spectral parameters
c
c  Original program had t as an input.  For C interface we need to 
c  center this on numpts to allow the arrays to be dynamically allocated
c  Hence t is dynamic.
c      numpts=ifix(t/dt+1.5)
      nft=npowr2(numpts)
c  because we give numpts as an input we must require it to be a 
c  power of 2 in this procedure.  Hence we have to trap this condition
      if(nft.ne.numpts) then
            ierr=-2
            return
      endif
      nfpts=nft/2+1
      fny=1./(2.*dt)
      delf=2.*fny/float(nft)
      t=dt*real(nft)
c
c     set up some computational parameters
c          specifying the type of response
c          requested.
c
      p = cmplx(pr,0.)
c original code used this alphanumeric character
c changed by glp to int as simpler to deal with in 
c C++ interface to this procedure
c      if ( complt(1:1) .eq. 'f' )  then
c this is the replacement
      if(fullrsp.ne.0) then
         rvb = allrvb
       else
         rvb = norvb
       endif
c-- original code
c      if ( modcnv(1:1) .eq. 'n' ) then
c       cnv = prmphs
c       else
c       cnv = allphs
c      endif
c--conversion frozen to answer yes
      cnv=allphs
c
c
c
c     compute q, alfa, and beta at 1 hz for absorbtion band
c
      t1 = 1.0d04
      wq = twopi
      do 5 i = 1, nlyrs
         qa = qpm(i)
         qb = qsm(i)
         t2 = ta(i)
         alfa(i) = alfm(i) * vabm(wq,t1,t2,qa)
         t2 = tb(i)
         beta(i) = betm(i) * vabm(wq,t1,t2,qb)
         qa = qabm(wq,t1,t2,qa)
         qb = qabm(wq,t1,t2,qb)
         alfa(i) = alfa(i)*( 1. + (0.,0.5)/qa)
         beta(i) = beta(i)*( 1. + (0.,0.5)/qb)
         cnvrsn(i) = cnv
         reverb(i) = rvb
         rho(i) = rhom(i)
 5       thik(i) = thikm(i)
      cnvrsn(0) = cnv
c  another reverb change for C - see above -
c      if ( complt(1:1) .ne. 'f' )  then
	   if(fullrsp.eq.0) then
         reverb(1) = onervb
       endif
c
      fr = cmplx(1.,0.)
      call ifmat(1,p,fr,nlyrs)
c
      do 10 i = 1, nfpts-1
         fr = cmplx(delf * ( i - 1 ), 0. )
         wq = twopi * fr
         do 6 j = 1, nlyrs
            qa = qpm(j)
            qb = qsm(j)
            t2 = ta(j)
            alfa(j) = alfm(j) * vabm(wq,t1,t2,qa)
            t2 = tb(j)
            beta(j) = betm(j) * vabm(wq,t1,t2,qb)
            qa = qabm(wq,t1,t2,qa)
            qb = qabm(wq,t1,t2,qb)
            alfa(j) = alfa(j)*( 1. + (0.,0.5)/qa)
            beta(j) = beta(j)*( 1. + (0.,0.5)/qb)
 6       continue
         call rcvrfn(p,fr,nlyrs,dvp,dvs,drp,drs,dts)
         u0(i) = dvp * (0.,-1.)*(-1.,0.)
         w0(i) = drp
         u1(i) = dvs
         w1(i) = drs * (0.,1.)
         tn(i) = dts
10    continue
      u0(nfpts) = (0.,0.)
      w0(nfpts) = (0.,0.)
      u1(nfpts) = (0.,0.)
      w0(nfpts) = (0.,0.)
      tn(nfpts) = (0.,0.)
      zpeak=0.0
      do 11 i=1,nfpts-1
            zamp=abs(u0(i))
            if(zamp.gt.zpeak) then
                  zpeak=zamp
            endif
   11 continue
      rwlev=zpeak*wlevel
      do 12 i=1,nfpts-1
            zamp=abs(u0(i))
            if(zamp.lt.rwlev) then
                  if(zamp.gt.0.0) then
                        rfr(i)=w0(i)*zamp/(rwlev*u0(i))
                  else
                        rfr(i)=w0(i)/cmplx(rwlev,rwlev)
                  endif
             else
                  rfr(i)=w0(i)/u0(i)
             endif
   12 continue
      rfr(nfpts) = (0.,0.)
c -- multiply by gaussian filter with width fsigma
c -- rf also requires a phase shift to translate to tshift instead of 0
      do 20 i=1,nfpts-1
            freal=delf * ( i - 1 )
            u0(i)=u0(i)*exp(-(freal*freal/(2.0*fsigma*fsigma)))
            u1(i)=u1(i)*exp(-(freal*freal/(2.0*fsigma*fsigma)))
            w0(i)=w0(i)*exp(-(freal*freal/(2.0*fsigma*fsigma)))
            w1(i)=w1(i)*exp(-(freal*freal/(2.0*fsigma*fsigma)))
            tn(i)=tn(i)*exp(-(freal*freal/(2.0*fsigma*fsigma)))
            rfr(i)=rfr(i)*exp(-(freal*freal/(2.0*fsigma*fsigma)))
            rfr(i)=rfr(i)*exp(cmplx(0.0,-twopi*freal*rftshift))
   20 continue
c compute inverse fourier transform of all 5 elements 
c in the work arrays.  Original code wrote sac files here
      call dfftr(u0,nft,-1,delf)
      call dfftr(w0,nft,-1,delf)
      call dfftr(u1,nft,-1,delf)
      call dfftr(w1,nft,-1,delf)
      call dfftr(tn,nft,-1,delf)
      call dfftr(rfr,nft,-1,delf)
      return
      end      
      
      subroutine rcvrfn(p,f,nlyrs,dvp,dvs,drp,drs,dts)
      integer nlyrs
      complex p,f
      complex dvp,dvs,drp,drs,dts
c
c        compute receiver function - free surface displacement from a
c        plane wave incident from below, on a stack of plane, parallel,
c        homogeneous layers
c        for a p, sv or sh wave incident
c        interface 0 is top of layer 1, a free surface,
c        layer n is half space
c        given frequency and phase slowness.
c
c          arguments...
c        psvsh = 1,2,3 for an incident p, sv or sh wave.
c  above comment confusing.  I think it is this.  dvp and dvr
c  are amplitude and phase for f for incident P wave on vertical (z)
c  and radial (r).  dvs and drs are rfor incident s.   dts is 
c  SH term
c
c        f,p - prescribed freq (hz) & horizontal phase slowness (c is
c            not restricted to be greater than alfa or beta)
c            both may be complex
c
c        passed in common /model/
c        alfa,beta,qp,qs,rho and thik contain the medium properties for
c
c        nlyrs - total number of layers, layer nlyrs is
c            the half space
c
c
c
c        commons and declarations
c
c
      include 'kennet.inc'   
c
c        complex declarations
c
      complex i,zero,one,w
      complex t11,t12,t21,t22,l11,l12,l21,l22,tsh,lsh
      complex*16 det
      complex x11,x12,x21,x22,y11,y12,y21,y22,xsh,ysh
      complex tnupp,tnups,tnusp,tnuss,tnush
      complex rndpp,rndps,rndsp,rndss,rndsh
      complex phtp,phts,phtpp,phtps,phtss
      real twopi
      integer lyr,nif,cnvnif
      external cphs
      complex cphs
      data twopi/6.2831853/,i,zero/(0.,1.),(0.,0.)/,one/(1.,0.)/
c
c
      w = twopi*f
c
c     handle the special case of a half space
c
      if ( nlyrs .eq. 1 ) then
         dvp = dvpfs
         dvs = dvsfs
         drp = drpfs
         drs = drsfs
         dts = dtshfs
         return
       endif
c
c        initialize tup and rdown matricies for the stack with
c        bottom interface matricies
c
      nif = nlyrs-1
      cnvnif = cnvrsn(nif)
      if ( cnvnif .eq. allphs ) then
         tnupp = tupp(nif)
         tnuss = tuss(nif)
         tnups = tups(nif)
         tnusp = tusp(nif)
         tnush = tush(nif)
         rndpp = rdpp(nif)
         rndss = rdss(nif)
         rndps = rdps(nif)
         rndsp = rdsp(nif)
         rndsh = rdsh(nif)
       else if ( cnvnif .eq. prmphs ) then
         tnupp = tupp(nif)
         tnuss = tuss(nif)
         tnups = zero
         tnusp = zero
         tnush = tush(nif)
         rndpp = rdpp(nif)
         rndss = rdss(nif)
         rndps = zero
         rndsp = zero
         rndsh = rdsh(nif)
       else if ( cnvnif .eq. cnvphs ) then
         tnups = tups(nif)
         tnusp = tusp(nif)
         tnupp = zero
         tnuss = zero
         tnush = tush(nif)
         rndps = rdps(nif)
         rndsp = rdsp(nif)
         rndpp = zero
         rndss = zero
         rndsh = rdsh(nif)
       endif
c
c        now do the  bottom up recursion for tup and rdown
c
      do 10 lyr = nlyrs-1, 2, -1
         nif = lyr - 1
c
c        use the two way phase delay through the layer
c        to/from the next interface
c
         phtp = cphs( -i*w*xi(lyr)*thik(lyr) )
         phts = cphs( -i*w*eta(lyr)*thik(lyr) )
         phtpp = phtp * phtp
         phtps = phtp * phts
         phtss = phts * phts
         rndpp = rndpp * phtpp
         rndss = rndss * phtss
         rndps = rndps * phtps
         rndsp = rndsp * phtps
         rndsh = rndsh * phtss
         tnupp = tnupp * phtp
         tnuss = tnuss * phts
         tnups = tnups * phtp
         tnusp = tnusp * phts
         tnush = tnush * phts
	 stnupp(lyr) = tnupp
	 stnups(lyr) = tnups
	 stnusp(lyr) = tnusp
	 stnuss(lyr) = tnuss
	 stnush(lyr) = tnush
	 srndpp(lyr) = rndpp
	 srndps(lyr) = rndps
	 srndsp(lyr) = rndsp
	 srndss(lyr) = rndss
	 srndsh(lyr) = rndsh
c
c        form the reverberation operator for the layer
c
         cnvnif = cnvrsn(nif)
         if ( cnvnif .eq. allphs ) then
            t11 = rupp(nif)
            t22 = russ(nif)
            t12 = rups(nif)
            t21 = rusp(nif)
            tsh = rush(nif)
          else if ( cnvnif .eq. prmphs ) then
            t11 = rupp(nif)
            t22 = russ(nif)
            t12 = zero
            t21 = zero
            tsh = rush(nif)
          else if ( cnvnif .eq. cnvphs ) then
            t12 = rups(nif)
            t21 = rusp(nif)
            t11 = zero
            t22 = zero
            tsh = rush(nif)
          endif
         if ( reverb(lyr) .eq. allrvb ) then
            l11 = one - (rndpp*t11 + rndps*t21)
            l22 = one - (rndsp*t12 + rndss*t22)
            l12 = - (rndpp*t12 + rndps*t22)
            l21 = - (rndsp*t11 + rndss*t21)
            det = ( l11*l22 - l12*l21 )
            l12 = -l12/det
            l21 = -l21/det
            t11 = l11/det
            l11 = l22/det
            l22 = t11
            lsh = one / ( one - rndsh*tsh )
         else if ( reverb(lyr) .eq. onervb ) then
            l11 = one + (rndpp*t11 + rndps*t21)
            l22 = one + (rndsp*t12 + rndss*t22)
            l12 =  (rndpp*t12 + rndps*t22)
            l21 =  (rndsp*t11 + rndss*t21)
            lsh = one + rndsh*tsh
         else if ( reverb(lyr) .eq. norvb ) then
            l11 = one
            l22 = one
            l12 = zero
            l21 = zero
            lsh = one
          endif
c
c        now finish the recursion, adding the next interface
c
         if ( cnvnif .eq. allphs ) then
            x11 = tupp(nif)
            x22 = tuss(nif)
            x12 = tups(nif)
            x21 = tusp(nif)
            xsh = tush(nif)
            y11 = rdpp(nif)
            y22 = rdss(nif)
            y12 = rdps(nif)
            y21 = rdsp(nif)
            ysh = rdsh(nif)
          else if ( cnvnif .eq. prmphs ) then
            x11 = tupp(nif)
            x22 = tuss(nif)
            x12 = zero
            x21 = zero
            xsh = tush(nif)
            y11 = rdpp(nif)
            y22 = rdss(nif)
            y12 = zero
            y21 = zero
            ysh = rdsh(nif)
          else if ( cnvnif .eq. cnvphs ) then
            x12 = tups(nif)
            x21 = tusp(nif)
            x11 = zero
            x22 = zero
            xsh = tush(nif)
            y12 = rdps(nif)
            y21 = rdsp(nif)
            y11 = zero
            y22 = zero
            ysh = rdsh(nif)
          endif
c
         t11 = l11*tnupp + l12*tnusp
         t22 = l21*tnups + l22*tnuss
         t21 = l21*tnupp + l22*tnusp
         t12 = l11*tnups + l12*tnuss
         tsh = lsh * tnush
c
c        tnupp = tupp(nif)*t11 + tups(nif)*t21
c        tnuss = tusp(nif)*t12 + tuss(nif)*t22
c        tnups = tupp(nif)*t12 + tups(nif)*t22
c        tnusp = tusp(nif)*t11 + tuss(nif)*t21
         tnupp = x11*t11 + x12*t21
         tnuss = x21*t12 + x22*t22
         tnups = x11*t12 + x12*t22
         tnusp = x21*t11 + x22*t21
         tnush = xsh * tsh
c
c        t11 = l11*tdpp(nif) + l21*tdsp(nif)
c        t12 = l11*tdps(nif) + l21*tdss(nif)
c        t21 = l12*tdpp(nif) + l22*tdsp(nif)
c        t22 = l12*tdps(nif) + l22*tdss(nif)
         t11 = l11*x11 + l21*x12
         t12 = l11*x21 + l21*x22
         t21 = l12*x11 + l22*x12
         t22 = l12*x21 + l22*x22
         tsh = lsh * xsh
         l11 = rndpp*t11 + rndps*t21
         l12 = rndpp*t12 + rndps*t22
         l21 = rndsp*t11 + rndss*t21
         l22 = rndsp*t12 + rndss*t22
         lsh = rndsh * tsh
c        rndpp = rdpp(nif) + tupp(nif)*l11 + tups(nif)*l21
c        rndss = rdss(nif) + tusp(nif)*l12 + tuss(nif)*l22
c        rndps = rdps(nif) + tupp(nif)*l12 + tups(nif)*l22
c        rndsp = rdsp(nif) + tusp(nif)*l11 + tuss(nif)*l21
         rndpp = y11 + x11*l11 + x12*l21
         rndss = y22 + x21*l12 + x22*l22
         rndps = y12 + x11*l12 + x12*l22
         rndsp = y21 + x21*l11 + x22*l21
         rndsh = ysh + xsh*lsh
c
10    continue
c
c        use the two way phase delay through the top layer
c
         phtp = cphs( -i*w*xi(lyr)*thik(lyr) )
         phts = cphs( -i*w*eta(lyr)*thik(lyr) )
         phtpp = phtp * phtp
         phtps = phtp * phts
         phtss = phts * phts
         tnupp = tnupp * phtp
         tnuss = tnuss * phts
         tnups = tnups * phtp
         tnusp = tnusp * phts
         tnush = tnush * phts
         rndpp = rndpp * phtpp
         rndss = rndss * phtss
         rndps = rndps * phtps
         rndsp = rndsp * phtps
         rndsh = rndsh * phtss
c
c        form the reverberation operator for the top layer
c
         cnvnif = cnvrsn(0)
         if ( cnvnif .eq. allphs ) then
            t11 = ruppfs
            t22 = russfs
            t12 = rupsfs
            t21 = ruspfs
            tsh = rushfs
          else if ( cnvnif .eq. prmphs ) then
            t11 = ruppfs
            t22 = russfs
            t12 = zero
            t21 = zero
            tsh = rushfs
          else if ( cnvnif .eq. cnvphs ) then
            t12 = rupsfs
            t21 = ruspfs
            t11 = zero
            t22 = zero
            tsh = rushfs
          endif
         if ( reverb(lyr) .eq. allrvb ) then
            l11 = one - (rndpp*t11 + rndps*t21)
            l22 = one - (rndsp*t12 + rndss*t22)
            l12 = - (rndpp*t12 + rndps*t22)
            l21 = - (rndsp*t11 + rndss*t21)
            det = ( l11*l22 - l12*l21 )
            l12 = -l12/det
            l21 = -l21/det
            t11 = l11/det
            l11 = l22/det
            l22 = t11
            lsh = one / ( one - rndsh*tsh )
         else if ( reverb(lyr) .eq. onervb ) then
            l11 = one + (rndpp*t11 + rndps*t21)
            l22 = one + (rndsp*t12 + rndss*t22)
            l12 =  (rndpp*t12 + rndps*t22)
            l21 =  (rndsp*t11 + rndss*t21)
            lsh = one + rndsh*tsh
         else if ( reverb(lyr) .eq. norvb ) then
            l11 = one
            l22 = one
            l12 = zero
            l21 = zero
            lsh = one
          endif
c
c        now add the free surface displacement
c
         t11 = l11*tnupp + l12*tnusp
         t22 = l21*tnups + l22*tnuss
         t21 = l21*tnupp + l22*tnusp
         t12 = l11*tnups + l12*tnuss
         tsh = lsh*tnush
         dvp = dvpfs*t11 + dvsfs*t21
         dvs = dvpfs*t12 + dvsfs*t22
         drp = drpfs*t11 + drsfs*t21
         drs = drpfs*t12 + drsfs*t22
         dts = dtshfs*tsh
c
c
c
      return
      end
