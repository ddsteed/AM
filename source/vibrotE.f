C     Program  vibrotE
c 
c     This program is to calculate the vibrational/ro-vibrational
c   energies, or scattering threshold energies and channel energies.
c 
c     Written  by  Dr. Weiguo Sun  on  05/01/1997 
c     Modified by  Dr. Weiguo Sun  and  Shilin Hou  
c                                  in  05/02/2000 ---  07/12/2000
c     Modified by  Dr. Weiguo Sun  on  05/05/2001 ---  08/05/2001
c----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      dimension E(200),Eryd(200),Ethre(200),Echan(200),Eev(200)
      dimension Evr(200,200),Eau(200),Evb(200),av(200)
      dimension nj(200)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectrb/ Be,alphae,gamae,Der,betae,bryd,wvib0(40)
      common /spectr0/ w0,we0,wei,wex,wey,wez,wet,wes,wer,v0,r0
      common /spectr1/ Bee,ale,gae,eta3,eta4,eta5,eta6,eta7
      common /spectr2/ Dee,bete,xsi2,xsi3,xsi4,xsi5,xsi6,xsi7
      common /Eswitch/ aye,aze,ate,ase,are
      common /Eswitc1/ abe,aae,age,ae3,ae4,ae5,ae6,ae7
      common /Eswitc2/ ade,abt,ax2,ax3,ax4,ax5,ax6,ax7
      common /Eswitc3/ bxe,bye,bze,bte,bse,bre
      common /Eswitc4/ cxe,cye,cze,cte,cse,cre
      common /bdt/ ak,beta,br4,br6,br8,br10,br12,br14
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /LA1/ Deinp,Ner,mset,mcst,mtd,mdp,mv,nvs(20,80)
      common /fncom/ xmean(30),gf(30)
c----------------------------------------------------------------------
c
c On input:
c
c  mtyp = 1, use input vib-rot constants to cal. Ev & Evj;
c       = 2, use INPUT Ev & d_V(R)/d_R to cal. fn, {Ev} & Evj.
c       = 3, calc. DISSOCIATION energy De from input Ev.
c
c   nw0 = 2, E(v) = w0 + (We + We0)*(v+0.5) - WeXe*(v+0.5)**2 + ...
c       = 1, E(v) = w0 +  We'*(v+0.5) - WeXe*(v+0.5)**2 + ...
c       = 0, E(v) = We*(v+0.5) - WeXe*(v+0.5)**2 + ...
c   nvd = number of vib. states used in comparison with
c           the input Ev or energy differencies. nvd =< n_expt.
c           nvd is slightly smaller than n_expt which is the
c           no. of input Ev OR energy differencies for comparison.
c    kn = number of EXTRA vib. states used in comparison. 
c   imv = maximum vib. quantum number used to calculate De(inp)
c           using INPUT vibrational constants.
c   kin = 1, READ expt'l energies or their differencies;
c       = 0, do NOT read.
c  nlin = 1, cal. fn's from CALC. vib. energies using "broydn";
c       = 2, cal. fn's from INPUT vib. energies using "broydn";
c       = 0, do NOT cal. force constants fn's.
C--
c mset -- number of set of solution vectors  X;
C           usually, set  mset =< 20.
C            !  mset =<  mv - mcst  !
C mcst -- number of vib. constants in a vector X;
C           usually, set  mcst =< 7.    But as solve
C           for X using INPUT Ev's, 7 < mcst < mv.
C             A(mcst, mcst);  b(mcst, 1)
C   mv -- number of vib. states in A used to
C           solve linear equation    A * X = b
C             set  mv =< 95.
C
C meny = 0 :
C  Ev = We(v+1/2) - WeXe(v+1/2)**2 + WeYe(v+1/2)**3
C                 + WeZe(v+1/2)**4 + WeTe(v+1/2)**5 + ... + ...
C meny = 1 :
C  Ev = We(v+1/2) - WeXe(v+1/2)**2 + WeYe(v+1/2)**3
C                 - WeZe(v+1/2)**4 + WeTe(v+1/2)**5 - ... + ...
C
C  mtd -- Decide which VIB. states being used to solve A*X = E :
C       = 0, use VIB. states in successively increasing order;
C       > 0, use VIB. states particularly specified.
C
C  mdp -- The mdp_th set from (A*X=E) CONVERGES spectrum { Ev }.
C 
C         FOR  mtd = 0 : 
C nvs(i,1) -- value of INITITAL vib. quantum state
C             v+1 in the ith set of mset vectors, e.g.
C                 i  = 1, 2, 3, 4, 5   (mcst=5 states)
C           nvs(i,1) = 1, 3, 6, 8, 9
C              v_ini = 0, 2, 5, 7, 8
C          ! nvs(i,1) + mcst =< mv   ! 
C    
C    When mcst=6 & mtd = 0, the states used in a set are :
C  nvs(1,1)=1 :  v = 0, 1, 2, 3, 4, 5.
C  nvs(3,1)=6 :  v = 5, 6, 7, 8, 9,10.
C
C    For ith set, if mcst=7 & mtd > 0 , arbitrary states are used :
C  nvs(i,j)=2, 3, 6, 7, 14, 15.
C         v=1, 2, 5, 6, 13, 14.
C
C  gf(i) -- The initial guess for vibrational FORCE constants : 
C    [ They may be used to cal. gf(i) using 'broydn' ]
C      gf(i) < 0, i = 3, 5, ...;  gf(i) > 0, i = 4, 6, ...
C    
C       For nw0 = 2 :
C    X = ( W0, We+We0, WeXe, WeYe, WeZe, WeTe, WeSe, ...)
C             We from input;  We0 from calculation.
C       For nw0 = 1  ( We' = We+We0 ) :
C    X = ( W0, We', WeXe, WeYe, WeZe, WeTe, WeSe, ...)
C             We' = We+We0  is obtained from  A*X = E.
C       OR (for nw0 = 0) :
C    X = ( We, WeXe, WeYe, WeZe, WeTe, WeSe, WeRe )
C--
C    When check VIBrational energies using the vibrational
c  constants generated from the INPUT vib. energies by
c  linear algebraic method (A*X=B).
c
c  bxe, ..., bre = 0, Zero calc. constants :
c                    WeXi, WeYi, WeZi, WeTi, WeSi, WeRi;
c                = +1, Do NOT change sign of these constants;
c                = -1, Change the sign of these constants.
C--
c   Mvr = 1, for vibration; = 2 for rot-vibration; =0 for BOTH.
c    NE = The # of input SCATTERING energies
c    Nv = # of vibrational states considered
c                Nv=6 means that v=0,1,2,3,4,5
c  iniv = Quantum # of the initial vibrational reference state
c  jnir = Quantum # of the initial  rotational reference state
c                 (Usually take iniv=v0=0, jnir=j0=0)
c   Nw  = 0, input (WE,...,Be,...) in a.u.;
c       = 1, input (WE,...,Be,...) in cm-1, code will convert them.
c nj(i) = # of rotational states in the i_th vibrational state.
c   lda = the eigenvalue of the z component Lz of the electronic
c         angular momentum L whose z axis coinciding with the
c         molecular axis.  "lda" presents the coupling between the
c         electronic and nuclear motion of the diatomic system.
c
c  beta = WIDTH (adjustable) parameter of the potential.
c
c        E(i)  = The energies of the scattered particle (in eV)
c    We, WeXe  = Experimental vibrational energy constants (in Hartree)
c Be,alphae,Der= Experimental  rotational energy constants (in Hartree)
c
c   Using Sun's perturbation formulae to calculate energies Evj.
c
c----------------------------------------------------------------------
       rbohr = 0.529177249d0
       rydev = 13.60569809d0
 	  auev = 27.21139618d0
        aucm = 219474.6306d0
        amae = 1.6605655D-27/9.1093897D-31
c----------------------------------------------------------------------
c      amae=mass_unit/mass_e
C      amae=1.6605655D-27/9.1093897D-31=1822.9163
c----------------------------------------------------------------------
      read (5,*) mtyp
      read (5,*) nw0,nvd,kn,imv,kin,nlin
      read (5,*) mset,mcst,mv,meny,mtd,mdp
        if (mtd .eq. 0) then  
	    do j=1,1
	      read (5,*) ( nvs(i,j), i=1,mset )
	    enddo
        else 
	    do i=1,mset
	      read (5,*) ( nvs(i,j), j=1,mcst )
	    enddo
        endif
      read (5,*) ( gf(i), i=3,mset+2 )
      read (5,*) bxe, bye, bze, bte, bse, bre
      read (5,*) Mvr,NE,Nv,iniv,jnir,Nw
	   beta = 1.0
Cc    read (5,*) Re,De,beta
      read (5,*) Re,De
      read (5,*) ams,bms,am0
        if (am0 .eq. 0.0) am0 = ams*bms/(ams+bms)
        if (mtyp .eq. 1) nlin = 0 
          amu = am0 * amae
        if (Nw .gt. 0) then
          Re  = Re/rbohr
          De  = De/aucm
        endif
c
        Deinp = De
c
      read (5,*) (nj(k), k=1,Nv)
        write( 6,100) mtyp,Mvr,NE,Nv,iniv,jnir,Nw,beta
        write(35,100) mtyp,Mvr,NE,Nv,iniv,jnir,Nw,beta
        write(38,100) mtyp,Mvr,NE,Nv,iniv,jnir,Nw,beta
        write( 6,102) 
        write(35,102) 
        write(38,102) 
        write( 6,400) nw0,nvd,kn,imv,kin,nlin,Re,De,ams,bms,am0,amu
        write(35,400) nw0,nvd,kn,imv,kin,nlin,Re,De,ams,bms,am0,amu
        write(38,400) nw0,nvd,kn,imv,kin,nlin,Re,De,ams,bms,am0,amu
      do k=1,Nv
        write( 6,105) k-1, nj(k) 
        write(35,105) k-1, nj(k) 
        write(38,105) k-1, nj(k) 
      enddo
c
        read (5,*) lda
        write( 6,118) lda
        write(35,118) lda
        write(38,118) lda
c
      write( 6,122) mset,mcst,mv,meny,mtd,mdp,bxe,bye,bze,bte,bse,bre
      write(35,122) mset,mcst,mv,meny,mtd,mdp,bxe,bye,bze,bte,bse,bre
      write(38,122) mset,mcst,mv,meny,mtd,mdp,bxe,bye,bze,bte,bse,bre
C-
        if ( mtd .eq. 0 ) then
            write( 6,144) 
            write(35,144) 
            write(38,144) 
              ms0 = 1
          do j=1,ms0
            do k=1,mset
                write( 6,105) k, nvs(k,j)
                write(35,105) k, nvs(k,j)
                write(38,105) k, nvs(k,j)
                  nk0 = nvs(k,j) + mcst
              if (nk0 .gt. mv ) then
		    write( 6,132) k, nvs(k,j), mcst, mv
		    write(38,132) k, nvs(k,j), mcst, mv
		    write(0,*)  " Error in mcst !  CHECK mcst ! "
		      stop
              endif
            enddo
          enddo
C-
        else
            write( 6,146) 
            write(35,146) 
            write(38,146) 
            write(80,146) 
              ms0 = mset
          do k=1,ms0
            write( 6,155) k, ( nvs(k,j), j=1,mcst )
            write(35,155) k, ( nvs(k,j), j=1,mcst )
            write(38,155) k, ( nvs(k,j), j=1,mcst )
            write(80,155) k, ( nvs(k,j), j=1,mcst )
          enddo
        endif
c---
          write( 6,125)
          write(35,125)
          write(38,125)
        do k=3,mset+2
          write( 6,110) k, gf(k)
          write(35,110) k, gf(k)
          write(38,110) k, gf(k)
        enddo
c
          cxe = bxe
          cye = bye
          cze = bze
          cte = bte
          cse = bse
          cre = bre
c----------------------------------------------------------------------
          nj1 = nj(1)
        nj(1) = 1
c----------------------------------------------------------------------
        v0=DFLOAT(iniv)
        r0=DFLOAT(jnir)
      if (NE .gt. 0) then
        read (5,*) (E(i), i=1,NE)
          write(6,108) 
        do i=1,NE
          write(6,110) i, E(i) 
        enddo 
      endif
c----------------------------------------------------------------------
      read (5,*) We,WeXe,WeYe
      read (5,*) WeZe,WeTe,WeSe
      read (5,*) WeRe
      read (5,*) Be,alphae,gamae
      read (5,*) Der,betae
        write( 6,115) We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe,Be,alphae,
     #gamae,Der,betae
        write(38,115) We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe,Be,alphae,
     #gamae,Der,betae
c-
      read(5,*) aye,aze,ate,ase,are
      read(5,*) abe,aae,age,ae3,ae4,ae5,ae6,ae7
      read(5,*) ade,abt,ax2,ax3,ax4,ax5,ax6,ax7
c--- 
      if (Nw .gt. 0) then
         We  = We/aucm
        WeXe = WeXe/aucm
        WeYe = WeYe/aucm
        WeZe = WeZe/aucm
        WeTe = WeTe/aucm
        WeSe = WeSe/aucm
        WeRe = WeRe/aucm
          Be =   Be/aucm
       alphae= alphae/aucm
       gamae = gamae/aucm
         Der =  Der/aucm
       betae = betae/aucm
      endif
         Wee = We
c-
       gf(2) = We*We*amu
c----------------------------------------------------------------------
          w0 = 0.0
         we0 = 0.0
         wex = WeXe
         wey = WeYe
         wez = WeZe
         wet = WeTe
         wes = WeSe
         wer = WeRe
         Bee = Be
         ale = alphae
         gae = gamae
         Dee = Der
        bete = betae
c----------------------------------------------------------------------
         if (mtyp .eq. 2 .or. mtyp .eq. 3) goto 10
c---
C switch mtyp.eq.3 was added by Hou on 4th,July,2000
c----------------------------------------------------------------------
c  Cal. vibrational-rotational REFERENCE energies Eref for mtyp=1
c----------------------------------------------------------------------
        call  calEvj(0,1,nj,lda,200,Ethre,Evr)
          Erefv = Ethre(1)
          Erefr = Evr(1,1)
c----------------------------------------------------------------------
c  Calculate vib.-rotational energies Ev==Evb; Evj=Evr for mtyp=1
c----------------------------------------------------------------------
          nj(1) = nj1
        call  calEvj(0,Nv,nj,lda,200,Evb,Evr)
          goto 35
c======================================================================
c  Cal. force constants and vib. constants used in Ev/Evj for mtyp=2
c----------------------------------------------------------------------
   10   call  VRener(Nv)
c----------------------------------------------------------------------
c     FOR  De_inp > 0.0 :
c  "VRener"[-> 'guessFn'-> 'vjcoef'-> 'convj'] -- > "calEvj" --> 
c  --> "getEvDe"(read experimental Ev's) --> ...
c
c----------------------------------------------------------------------
c     FOR  De_inp = 0.0 :
c  "getEvDe"(read experimental Ev's) --> "calVIBconst" -> "linersolv"
c  (solve A*X=E) --> "vibev"(cal. Ev) --> "guessFn"(cal. initial guess
c  of fn's) --> "broydn"(solve for fn's) --> "vibconst"(VIBrational
c   constants) --> "vibev"  &  "FindnumDe"(cal. De) .
c----------------------------------------------------------------------
            IF (mtyp .eq. 3 .and. Deinp .eq. 0.0)  THEN
              call  getEvDe(Nv,0,lda,200,Evb)
                goto  900
            ENDIF
c----------------------------------------------------------------------
c  Cal. vibrational-rotational REFERENCE energies  Eref for mtyp=2
c----------------------------------------------------------------------
        call  calEvj(0,1,nj,lda,200,Ethre,Evr)
          Erefv = Ethre(1)
          Erefr = Evr(1,1)
c----------------------------------------------------------------------
c  Calculate vib.-rotational energies Ev==Evb; Evj=Evr  for mtyp>=2
c----------------------------------------------------------------------
   20     nj(1) = nj1
        call  calEvj(0,Nv,nj,lda,200,Evb,Evr)
c======================================================================
   35   if (Mvr .eq. 2) go to 50 
c----------------------------------------------------------------------
c  Calculate vibrational THREshold energies  Ethre.
c----------------------------------------------------------------------
          write(6,120) 
      do k=1,Nv 
	  kv=k-1   
	  av(k)=DFLOAT(kv) 
	  avn=av(k)
c--- 
c Ethre(k) are in Hartree
c        Ethre(k)=Wee*(avn-v0) - WeXe*(avn*(avn+1)-v0*(v0+1))    
c--- 
         Ethre(k)=Evb(k) - Erefv
c Eev(k) are in eV
	     Eev(k)=Ethre(k)*auev
c Convert energy from eV to Rydberg
          Eryd(k)=Eev(k)/rydev
c--- 
          dife=Eev(k) - Eev(k-1)
        write(6,130) kv,Eryd(k),Ethre(k),Eev(k),dife
      enddo
c
c----------------------------------------------------------------------
c  Calculate scattering (vibrational) CHANNEL energies  Echan.
c----------------------------------------------------------------------
      DO i=1,NE
c Convert energy from eV to Hartree (a.u.)
          Eau(i)=E(i)/auev
          write(6,140) Eau(i),E(i)
	  do k=1,Nv
	    Echan(k)=Eau(i) - Ethre(k)
            Eev(k)=Echan(k)*auev
          write(6,130) k-1,Echan(k),Eev(k)
        enddo
      ENDDO  
          write(6,*)
c----------------------------------------------------------------------
        if (Mvr .eq. 1) go to 900
C======================================================================
c  Calculate vibrational-rotational THREshold energies  Evr.
c----------------------------------------------------------------------
   50     write(6,150)
C======================================================================
           Evr0 = Erefr
      if (nj(1) .gt. 0) then
        do 65 k=1,Nv
             kv = k - 1
            njr = nj(k)
               Devr = 0.0
          do 60 i=1,njr
              ir = i - 1
             Evr(k,i) = Evr(k,i) - Evr0
           if (i .gt. 1) Devr = Evr(k,i) - Evr(k,i-1)
              Evrev = Evr(k,i)*auev
               Devr = Devr*auev
             write(6,160) kv,ir,Evr(k,i),Evrev,Devr
  60      continue
 	       write(6,*)
  65    continue
      endif
C----------------------------------------------------------------------
c  Calculate vibrational-rotational CHANNEL energies  Echan.
c----------------------------------------------------------------------
      IF (NE .gt. 0) THEN
        do 80 i=1,NE
	    Eau(i)=E(i)/auev
	    write(6,170) Eau(i),E(i)
	  do 75 k=1,Nv
	       kv = k-1
            njr = nj(k)
          do 70 j=1,njr
	         jr = j-1
	        Echan(k)=Eau(i) - Evr(k,j)
	          Eev(k)=Echan(k)*auev
            write(6,160) kv,jr,Echan(k),Eev(k)
  70      continue
            write(6,*) 
  75      continue
  80    continue
      ENDIF
C---
        write(6,200) 
          kk = 0
      do i=1,Nv
	    if ( Evb(i) .lt. Evb(i-1) .and. kk .eq. 0 ) then
            kk = 1
            write(6,*)
          endif
        write(6,210) i-1, Evb(i)
      enddo
C-
        write(6,220) 
          kk = 0
      do i=1,Nv
	    if ( Evb(i) .lt. Evb(i-1) .and. kk .eq. 0 ) then
            kk = 1
            write(6,*)
          endif
        write(6,230) i-1, Evb(i)*aucm
      enddo
C-
        write(6,240) 
      do i=1,Nv
         dif = Evb(i) - Evb(i-1)
           if (i .eq. 1) dif = 0.0
        write(6,230) i-1, dif*aucm
      enddo
C----------------------------------------------------------------------
  100 format(///16x,'**=  Output of "vibrotE.f"  =** ',/
     #/9x,'Switches to calculate VIB-ROTational energies :',/
     #/6x,'[ mtyp = 1, use input vib-rot consts to cal. Ev & Evj;',
     #/13x,'= 2, use INPUT Ev & d_V(R)/d_R to cal. fn, {Ev} & Evj.',
     #/6x,'       = 3, calc. DISSOCIATION energy De from input Ev.',/
     #/6x,'            ! Run both mtyp=1, 2 to compare Ev(max) ! ',/
     #/6x,'   Mvr = 1, for vib.; = 2 for rot-vib.; = 0, for BOTH.',
     #/6x,'    NE = The number of input SCATTERING energies.',
     #/6x,'    Nv = The number of vibrational states considered.',
     #/6x,'              Nv=6 means that v=0,1,2,3,4,5 ',
     #/6x,'  iniv = Quantum no. of vibrational reference state.',
     #/6x,'  jnir = Quantum no. of  rotational reference state.',
     #/6x,'    Nw = 0, input (WE,...,Be,...) in a.u.; ',
     #/6x,'       = 1, input in cm-1, code will convert them.     ',
     #/6x,'  beta = WIDTH parameter of a Rydberg-like potential. ]',
     #//10x,'mtyp  Mvr  NE  Nv  iniv  jnir  Nw   beta ',
     # /10x,i3,3x,i2,3x,i2,2x,i2,2x,i3,3x,i3,3x,i2,1x,f10.7//,)
 102  format(6x,'[ nw0 = 2, E(v) = w0 + (We + we0)*(v+0.5) - ... ',
     #/20x,'We from INPUT;  w0 & we0 from calculations. ', 
     #/6x,"      = 1, E(v) = w0 + We'*(v+0.5) - WeXe*(v+0.5)**2 +...",
     #/20x,"We' = We_solv + we0   from solving   A*X = E . ", 
     #/6x,'      = 0, E(v) = We*(v+0.5) - WeXe*(v+0.5)**2 + ... ',
     #/6x,'  nvd = number of vib. states used in comparison with ',
     #/14x,'the input Ev or energy differencies. nvd =< n_expt. ',
     #/14x,'nvd may be slightly smaller than n_expt which is the ',
     #/14x,'no. of input Ev OR energy differencies for comparison. ',
     #//14x,'   The mv below MUST be  mv =< nvd ;    nvd =< 95.   ',
     #//7x,'  kn = number of EXTRA vib. states used in comparison.',
     #/6x,'  imv = maximum vib. quantum number used to calculate ',
     #/14x,'De(inp) using INPUT vibrational constants. ',
     #/6x,"  kin = 1, READ expt'l energies or their differencies; ",
     #/6x,'      = 0, do NOT read.',
     #/6x," nlin = 1, cal. fn's from CALC. vib. Ev using 'broydn';",
     #/6x,"      = 2, cal. fn's from INPUT vib. Ev using 'broydn';",
     #/6x,"      = 0, do NOT cal. force constants fn's.  ",
     #/6x,"             Set  nlin = 0  for mtyp = 1  ! ",19x,']', )
 400  format(//20x,'nw0  nvd   kn    imv  kin  nlin',
     #/20x,i2,3x,i3,2x,i3,3x,i4,2x,i2,4x,i2//
     #/6x,'{   Re = Equilibrium internuclear distance in ao. ',
     #/6x,'    De = The dissociation energy of AB in a.u.    ',/
     #/6x,'           Input De = 0  TO GUESS initial De !  ',/
     #/6x,'  Amas = The MASS of atom A in amu.               ',
     #/6x,'  Bmas = The MASS of atom B in amu.               } ',
     #//9x,'                                    Reduced_mass_of_AB ',
     # /3x,'Re(ao)',3x,'De(a.u)     Amas       Bmas',
     #'     (in amu)    (in a.u) ',
     # /1f9.5,1f10.6,1x,1f10.6,1x,1f10.6,1x,1f10.6,1f14.6,
     #//3x,'   nj = is no. of rotational states in each vib. state :',
     #//23x,'  v      nj  ',/)
  420 format(5x,10(i2,3x) )
  105 format(23x,i3,4x,i4)
  108 format(//15x,'The input SCATTERING energies : ',
     #//23x,'  i      E(i;eV) ',/)
  110 format(23x,i3,3x,f10.5)
  115 format(//2x,'   The input ro-vibrational constants are : '//,
     #7x,'  We ',11x,'WeXe',11x,' WeYe',/3(1PE16.8),//
     #7x,'WeZe ',11x,'WeTe',11x,' WeSe',/3(1PE16.8),//
     #7x,'WeRe ',/1PE16.8,//
     #7x,'  Be ',10x,'Alphae',10x,'Gammae',/3(1PE16.8),//
     #7x,' Der ',11x,'betae',/2(1PE16.8),/)
  118 format(//5x,'The eigenvalue "lda" of the z component Lz of the',
     #/5x,'electronic angular momentum L whose z axis coinciding',
     #/5x,'with the molecular axis.  "lda" presents the coupling',
     #/5x,'between the electronic and nuclear motion of the molecule',
     #//5x,'  lda is used in  [ J*(J+1) - lda*lda ]  for ROTation ',
     #//5x,'                    lda =',i3)
  120 format(//15x,'** Vibrational threshold & channel energies **'//,
     #18x,'The vibrational threshold energies are : ',//,
     #3x,'v     Ethre(Rydberg)    Ethre(a.u.)      Ethre(eV) ',
     #'   Dif(v-1,v; eV)',/)
  122 format(///6x,'[ mset = number of set of solution vectors  X; ',
     #/17x,'usually, set  mset =< 20. ',
     #/17x,'!  mset =<  mv - mcst  ! ',
     #//6x,'  mcst = NUMBER of VIB. constants in a vector X; ',
     #/15x,'usually, set  mcst =< 7.      BUT as solve ',
     #/15x,"for X using INPUT Ev's,  7 =< mcst < mv. ",
     #/21x,'A(mcst, mcst);  b(mcst, 1) ',
     #//18x,'Suggest  mcst = 7  for  nw0 = 0 !  ',
     #//6x,'    mv = number of vib. states in A used to ',
     #/17x,'solve linear equation    A * X = b ',/
     #/19x,'set  mv =< nvd ;    nvd =< 95.   ',/
     #/10x,'mset + mcst =< mv;  otherwise calculation DIVERGES. ',/
     #/12x,'As  meny = 0 : ',
     #/8x,'Ev = We(v+1/2) - WeXe(v+1/2)**2 + WeYe(v+1/2)**3 ',
     #/18x,'     + WeZe(v+1/2)**4 + WeTe(v+1/2)**5 + ... + ... ',
     #/12x,'As  meny = 1 : ',
     #/8x,'Ev = We(v+1/2) - WeXe(v+1/2)**2 + WeYe(v+1/2)**3 ',
     #/18x,'     - WeZe(v+1/2)**4 + WeTe(v+1/2)**5 - ... + ... ',
     #//11x,'   For nw0 = 2 : ',
     #/11x,'X = (W0,We+We0,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe, ...)  ',
     #/11x,'          We from INPUT;  We0 from calculations. ',
     #/11x,'   For nw0 = 1 : ',
     #/11x,"X = (W0,We',WeXe,WeYe,WeZe,WeTe,WeSe,WeRe, ...)  ",
     #/11x,"          We'= We+We0  from solving  A*X = E . ",
     #/11x,'   OR for nw0 = 0 : ',
     #/11x,'X = ( We, WeXe, WeYe, WeZe, WeTe, WeSe, WeRe, ...)  ',
     #//8x,'mtd -- Decide which VIB. states being used to solve ',
     #'A*X = E : ',/13x,'= 0, use VIB. states ',
     #'in successively increasing order; ',
     #/13x,'> 0, use VIB. states particularly specified.',13x,']',
     #//11x,'mdp -- The mdp_th set from (A*X=E) CONVERGES ',
     #'spectrum { Ev }',//21x,'  mset + mcst =< mv',/
     #/21x,'mset  mcst   mv  meny  mtd  mdp ',
     #/21x,i3,3x,i3,3x,i3,3x,i2,3x,i2,3x,i2/
     #//4x,'[ When check VIBrational energies using the vibrational',
     # /4x,'  constants generated from the INPUT vib. energies by ',
     # /4x,'  linear algebraic method (A*X=B). ',
     #//4x,'  bxe, ..., bre = 0, Zero calc. constants : ',
     # /4x,'                     WeXi,WeYi,WeZi,WeTi,WeSi,WeRi; ',
     #/20x,'= +1, Do NOT change sign of these constants; ',
     #/20x,'= -1, Change the sign of these constants.    ] ',
     #//10x,' bxe   bye   bze   bte   bse   bre ',/10x,6(f4.1,2x),
     #///4x,'[ For mtd = 0 :    ',
     #/9x,'nvs(i,1) = value of initial vib. quantum state ',
     #/20x,'v+1 in the i_th set of mset vectors, e.g. ',
     #/15x,'       i = 1, 2, 3, 4, 5   (if mset=5 sets) ',
     #/15x,'nvs(i,1) = 1, 3, 6, 8, 9   ',
     #/15x,'   v_ini = 0, 2, 5, 7, 8   ',/
     #/16x,' !  nvs(i,1) + mcst =< mv  !',
     #///6x,'For mtd > 0 : ',
     #/9x,'Using mcst arbitrary VIB. states in i_th set of mset sets',
     #/17x,'nvs(i,j) = 2, 3, 6, 7, 14, 15.  ',
     #/17x,'       v = 1, 2, 5, 6, 13, 14.  ',19x,']',/)
  125 format(///7x,'Initial GUESSES of vibrational FORCE constants :',
     #//6x,"      Be SURE you have  mset  INPUT  f_n's  !      ",
     #//6x,"(Used to cal. f_n's using 'broydn' & We for ksh=0) ",
     #/18x," f_n(i) < 0, i = 3, 5, ... ; ",
     #/18x," f_n(i) > 0, i = 4, 6, ... . ",
     #//23x,'  i      f_n(i) ',/)
  130 format(x,i3,4f16.7)
  132 format(//3x," k  nvs(k,j)   mcst   mv ",
     #/3x,i2,4x,i2,5x,i3,4x,i2,
     #//3x,'    nvs(k,j) + mcst > mv  !    CHANGE mcst !  ',/)
  140 format(///1x,'For scattering energy =',f12.6,'  au  =',f12.6,
     #' eV',//' The channel energies (Ev = Kv**2 = 2*E) are : ',
     #//3x,'v      Echan(a.u.)      Echan(eV)  ',/)
  144 format(//23x,'    For  mtd = 0 :',//23x,'  i    nvs(i,1) ',/)
  146 format(//6x,'   For  mtd > 0 :  ',
     #/6x,' i         nvs(i,j)  j = 1, mcst ',
     #/6x,3(1H-),3x,32(1H-),)
C
  150 format(/6x,' * Ro-vibrational threshold & channel energies *'//,
     #9x,' The ro-vibrational threshold energies are : ',//,
     #3x,'      v   j   Ethre(a.u.)   Ethre(eV)   Delta_Evj(eV) ',/)
  155 format(/6x,i2,4x,10(i2,x),)
  160 format(7x,i3,i4,3f13.7)
  170 format(///2x,' For scattering energy =',f12.6,'  au  =',f12.6,
     #' eV',//,6x,' The channel energies (Evj = Kvj**2 = 2*E) are :',
     #//3x,'      v   j   Echan(a.u.)   Echan(eV)  ',/)
  200 format(//2x,'* Vibrational energies *'//,
     #3x,'  v      Evib(a.u.) ',/)
  210 format(3x,i3,1PE18.8)
  220 format(//2x,'* Vibrational energies *'//,
     #3x,'  v       Evib(cm-1) ',/)
  230 format(3x,i3,F18.6,F16.6)
  240 format(//3x,'* Vibrational energies *'//,
     #3x,'  v    Delta_Evib(cm-1) ',/)
C 240 format(//12x,'* Vibrational energies *'//,
C    #3x,'  v       Evib(cm-1)    Delta_Evib(cm-1) ',/)
C----------------------------------------------------------------------
  900   stop
      end
C
C==========
      subroutine  VRener(nv)
c-----------------------------------------------------------------
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /Eswitch/ aye,aze,ate,ase,are
      common /Eswitc1/ abe,aae,age,ae3,ae4,ae5,ae6,ae7
      common /Eswitc2/ ade,abt,ax2,ax3,ax4,ax5,ax6,ax7
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectrb/ Be,alphae,gamae,Der,betae,bryd,wvib0(40)
      common /spectr0/ w0,we0,wei,wex,wey,wez,wet,wes,wer,v0,r0
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /bdt/ ak,beta,br4,br6,br8,br10,br12,br14
      common /LA1/ Deinp,Ner,mset,mcst,mtd,mdp,mv,nvs(20,80)
      common /hnk/ hs2,Nryd,ksh
      common /forcecon/ f1,f2,f3,f4,f5,f6,f7,f8
      common /forceco1/ f9,f10,f11,f12,f13,f14,f15
      common /fncom/ xmean(30),gf(30)
      common /fnco1/ gg1(300),gg2(300),ge1(300),ge2(300)
      common /fnco2/ gg(300),ge(300),ggf(300)
c-----------------------------------------------------------------
c    Read data :
c
c     Use following V(R) to calculate n_th force constants f_n's.
c      |
c      | = 1, use Vmorse as V(R).
c      |
c      | = 2, use  Vryd(R) = -De *[ 1 + a*x ] *exp(-a*x)         1.)
c      |                a = b * dsqrt( amu/De ) ;   x = R - Re.
c      |                b = We + bryd ;  originally, b=We.
c      |--------------------------------------------------------
c      | The GENERALIZED Rydberg potentials are defined as :
c      |
c      | = 3, use
c      |      func(R)=[exp(-a*x)/beta - (a+bryd)*x]*exp(-a*x)
c      |      Vryd(R)= - De*beta*func(R)                         2.)
c      |           a = We * dsqrt( amu/De )
c      |
c      |------------------
c      | = 4, use
c      |      Uryd(R)=-De*[1.0 + (a-bryd)*x]*exp(-a*x)           3.)
c      |         a = bryd + dsqrt(bryd*bryd + f2a/De)
c      |
c      |------------------
c      | = 5, use
c      |      Uryd(R)=-De*[1.0 + a*x]*exp[-(a-bryd)*x]           4.)
c      |         a = dsqrt(bryd*bryd + f2a/De)
c      |------------------
c      | = 6, use
c      |      Umos(R)=De*[exp(-2*B*x) - (2+bryd*x)*exp(-B*x)]    5.)
c      |      Vcomb(R) = 0.5*[ Uryd(R) + Umos(R) ]
c      |             B = We * dsqrt{ amu/(2*De) }
c      |
c      |------------------
c      | = 7, use
c      |      Vmos(R)= De*[exp(-2*B*x) - 2*exp(-B*x)]
c      |      Vryd(R)=-De*[1.0 + a*x]*exp(-a*x)
c      |      V(R)=Vryd(R) + b0*[ Vryd(R) - Vmos(R) ]            6.)
c      |-----
c      | = 8, use
c      |      V(R)=Vmos(R) + b0*[ Vmos(R) - Vryd(R) ]            7.)
c      |        where  b0=bryd**Mryd  or  b0=(bryd/Re)**Mryd
c      |          and [...] has the same sign as the 1st term.
c Ny = |------------------
c      | = 9, use
c      |      Va(R)=Vcomb(R) + b0*[ Vryd(R) - Vmos(R) ]          8.)
c      |
c      |------
c      | = 10, use
c      |      Vb(R)=Vcomb(R) + b0*[ Vmos(R) - Vryd(R) ]          9.)
c      |
c      |------------------
c      | = 11, use
c      |      U(R) = [ Um(R) + Ur(R) ]/2.0                      10.)
c      |  where
c      |    Um(R) = De*[exp(-2*bt*x) - (2+bt*x)*exp(-bt*x)]
c      |    Ur(R) = -De*[1 + (af-bt)*x]*exp(-af*x)
c      |    bt = 0.5*(-cb + dsqrt(cb**2 + 2.0*ff2/De) )
c      |    cb = 0.5*( 1.0-sqrt(2.0) )*bryd
c      |
c      |------
c      | = 12, use
c      |      U(R) = [ Um(R) + Ur(R) ]/2.0                      11.)
c      |  where
c      |    Um(R) = De*[exp(-2*bt*x) - (2+bt*x)*exp(-bt*x)]
c      |    Ur(R) = -De*[1 + (af-bt)*x]*exp(-af*x)
c      |    bt = 0.5*(-cb + dsqrt(cb**2 + 4.0*ff2*bk/De) )/bk
c      |    cb = 0.5*( 1.0 - ak)*bryd  ;
c      |    bk = 1.0 + ak*ak/2.0 ;        af = ak*bt .
c      |    ak = INPUT value which is close to sqrt(2.0)=1.4142
c      |
c      |------
c      | = 13, use
c      |      U(R) = [ Um(R) + Ur(R) ]/2.0                      12.)
c      |  where
c      |    Um(R) = De*[(1-b*x)*exp(-2*bt*x) - 2*exp(-bt*x)]
c      |    Ur(R) = -De*[1 + (af-b)*x]*exp(-af*x)
c      |    bt = 0.5*(-cb + dsqrt(cb**2 + 2.0*ff2/De) )
c      |    cb = ( 1.0 - sqrt(2.0)/2 )*b  ;     b == bryd .
c      |
c      |------
c      | = 18, use Eq.12) & PVM procedures to calculate the
c      |       INItial guesses of the FORCE constants fn's
c      |       by using VIBrational constants from AM: A*X = E.
c      |
c      |==================
c      | = 26, use Vpseudo_gaussian  as V(R).
c
c  Ner = number of force constants f's wanted.
c
c bryd - The variational constant used to adjust Vryd(R).
c Mryd = M = number of powers in Eq 7.) for Ny = 8 .
c      bryd, Mryd are meaningless for Vmorse & Vp_g.
c
c  hs -- Estimated initial stepsize used by code "dfridr".
c            Good range for  hs :  0.001 --> 4.0  for H_2.
c              hs1 for Vrydberg;   hs2 for V_p-g .
c
c   Mwe = 1, New  We = We + bryd ;  = 0, We = We.
c  Nryd = Number of f_n's calculated using Vryd(R) and
c           is used to cal. V(R->Rmax=inf).
c  Rmax = The maximum R value for  V(R->Rmax=inf).
c    ak = a value close to sqrt(2.0)=1.4142 for Ny=21.
c         ak is used in alpha = ak*beta  in potential : 
c           U(R) = [ Um(R) + Ur(R) ]/2.0
c Rless = The R value in    Rbig = Rmax +/- Rless.
c Dconv = The tolerance used to check if V(Rmax) = De.
c
c br4,br6,br8,br10,br12,br14 = The scaling switches for the
c long range force const f_n(long;R=Re) = d~n/(dr~n)[1/R**k]
c which will be added into f_n(R=Re) = V_rydberg~n(R=Re).
c   f_n(long;R=Re) = (-1.0)**n *(k+n-1)!/[(k-1)!*R**(k+n)]
c which are important for Wande-Waals molecules and
c quasi-stable molecules.
c    brn = 1.0, Add f_n(long;R=Re); = 1.0, do NOT add.
c
c ksh = 0, use the INPUT first guess for FORCE constants fn's ;
c     = 1, use fn's from Morse V(R);
c     = 2, use fn's from Rydberg-like V(R);
c     = 3, use fn's from P-G V(R).
c=========================================================================
c  Prepare CONVERSION factors :
c-------------------------------------------------------- 
         aotoA0 = 0.5291772490
         aJtoau = 0.229371d0
c-------------------------------------------------------- 
             read(5,*) Ner,Ny,bryd,Mryd,hs1,hs2
             read(5,*) Mwe,Nryd,Rmax,ak,Rless,Dconv
             read(5,*) br4,br6,br8,br10,br12,br14
             read(5,*) ksh
               ij = 6
	     do i=1,3
             write(ij,790) 
             write(ij,800) Mwe,Ner,Ny,ak,bryd,Mryd
             write(ij,810) ksh,hs1,hs2
               if (ij .eq. 38) write(38,812) Ny, bryd
             write(ij,820) aye,aze,ate,ase,are,abe,aae,age,
     # ae3,ae4,ae5,ae6,ae7,ade,abt,ax2,ax3,ax4,ax5,ax6,ax7,
     # br4,br6,br8,br10,br12,br14
               if (i .eq. 1) ij = 35
               if (i .eq. 2) ij = 38
	     enddo
c-------------------------------------------------------- 
C   For INPUT De=0, goto 400 to call "vjcoef" to 
C generate  coefficients "b1,b3,b4,b5,b6,b7,b8" which
C will be used in "FUNCV" later.
c-------------------------------------------------------- 
         IF ( Deinp .eq. 0.0 ) goto  400
c-------------------------------------------------------- 
c   Generate the data needed by V(Pseudo-Gaussian; R).
c     ff2 = amu*We*We, & betap=dsqrt(ff2/(2.0*De) )
c        are generated in  "Vpotdata".   
c-------------------------------------------------------- 
             call  Vpotdata
c-------------------------------------------------------- 
c  Calculate force constants ff == gg  using
c  the definitions :   f_n = [ d~nV(R) ]/[dR~n] 
c  and  the function code "dfridr" .
c-------------------------------------------------------- 
             call  guessFn(We)
c-------------------------------------------------------- 
                 alpha0 = We*dsqrt( amu/De )
                  beta0 = We*dsqrt( amu/(2.0*De) )
               write( 6,534) De, Re, alpha0, beta0
               write(35,534) De, Re, alpha0, beta0
c            do k=2,Ner+1
             do k=1,Ner+1
               write( 6,550) k, gg(k), ge(k)
               write(35,550) k, gg(k), ge(k)
             enddo
c---
               write( 6,536)
               write(33,540) De, We, bryd
               write(38,540) De, We, bryd
C            do k=1,Ner+1
             do k=1,20
                ery = ggf(k)/(aJtoau * aotoA0**k)
               write( 6,550) k, ggf(k), ery
               write(33,550) k, ggf(k), ery
               write(38,550) k, ggf(k), ery
             enddo
c--------------------------------------------------- 
c Print numerical "force constants"
c--------------------------------------------------- 
             write( 6,750) 
             write(35,750) 
             write(25,752) 
c          do k=2,Nryd
           do k=1,Nryd
               write( 6,550) k, gg1(k), gg(k),  gg2(k)
               write(35,550) k, gg1(k), gg(k),  gg2(k)
               write(25,550) k, gg1(k)/(aJtoau * aotoA0**k),
     #  gg(k)/(aJtoau * aotoA0**k),  gg2(k)/(aJtoau * aotoA0**k)
C               ggk =  gg(k)/(aJtoau * aotoA0**k)
             if (k .eq. Ner+1) then
               write( 6,*)
               write(35,*)
             endif
           enddo
             write( 6,760) 
             write(35,760) 
c          do k=2,Ner+1
           do k=1,Ner+1
             write( 6,550) k, ge1(k), ge(k), ge2(k)
             write(35,550) k, ge1(k), ge(k), ge2(k)
           enddo
c---------------------------------------------------
c  V(R->inf) = De = Sum_n  f_n*(R - Re)**n/n!
c---------------------------------------------------
                kk = 0
                kj = 0
                bk = 1.0
              vmsg = 0.0
              vmsl = 0.0
              vdeg = 0.0
              vdel = 0.0
                Rg = Rmax + Rless
                Rl = Rmax - Rless
            write( 6,730) De, Rmax, Rless, Dconv, Nryd
            write(35,730) De, Rmax, Rless, Dconv, Nryd
C         do k=2,Nryd
          do k=1,Nryd
              bk = bk*k
            vmg0 = vmsg
            vmsg = vmsg + gg1(k)*(Rg - Re)**k/bk
            vmsl = vmsl + gg1(k)*(Rl - Re)**k/bk
c
            vdg0 = vdeg
            vdeg = vdeg + ggf(k)*(Rg - Re)**k/bk
            vdel = vdel + ggf(k)*(Rl - Re)**k/bk
       write( 6,740) k, vmsg, vmsl, vdeg, vdel, vdeg-vdel
       write(35,740) k, vmsg, vmsl, vdeg, vdel, vdeg-vdel
       If ( k .gt.1 ) then
            if ( abs(vdeg - vdel) .lt. Dconv .and. kk .eq. 0 ) then
              write( 6,732)
              write(35,732)
                kk = 1
                kc = k
               vde = vdeg
            endif
       Endif
       If ( k .gt.1 ) then
            if ( abs(vdeg - vdg0) .lt. Dconv .and. kk .eq. 0 ) then
              write( 6,732)
              write(35,732)
                kk = 1
                kc = k
               vde = vdeg
            endif
       Endif
       If ( k .gt.1 ) then
            if ( abs(vmsg - vmsl) .lt. Dconv .and. kj .eq. 0 ) then
              write( 6,734)
              write(35,734)
                kj = 1
                km = k
               vdm = vmsg
            endif
       Endif

       If ( k .gt.1 ) then
            if ( abs(vmsg - vmg0) .lt. Dconv .and. kj .eq. 0 ) then
              write( 6,734)
              write(35,734)
                kj = 1
                km = k
               vdm = vmsg
            endif
       Endif
          enddo
            write( 6,742) kc, vde, km, vdm
            write(35,742) kc, vde, km, vdm
C=================================================
 400       IF (Ner .gt. 0) then
c-------------------------------------------------
c  Define NEW vibrational constant 'We' if wanted
c-------------------------------------------------
              f1 = gg(1)
              f2 = gg(2)
              f3 = gg(3)
              f4 = gg(4)
c
            if (Ner .gt. 3) f5 = gg(5)
            if (Ner .gt. 4) f6 = gg(6)
            if (Ner .gt. 5) f7 = gg(7)
            if (Ner .gt. 6) f8 = gg(8)
c
            if (Ner .gt. 7)  f9  = gg(9)
            if (Ner .gt. 8)  f10 = gg(10)
            if (Ner .gt. 9)  f11 = gg(11)
            if (Ner .gt. 10) f12 = gg(12)
            if (Ner .gt. 11) f13 = gg(13)
            if (Ner .gt. 12) f14 = gg(14)
            if (Ner .gt. 13) f15 = gg(15)
            if (Ner .gt. 14) f16 = gg(16)
            if (Ner .gt. 15) f17 = gg(17)
            if (Ner .gt. 16) f18 = gg(18)
            if (Ner .gt. 17) f19 = gg(19)
            if (Ner .gt. 18) f20 = gg(20)
c---
             if (Mwe .gt. 0) then
                 We = We + bryd
               write( 6,745) Wee, We
             endif
c-------------------------------------------------
c  Calculate coefficients for ro-vib constants
c-------------------------------------------------
             call  vjcoef(We)
c-------------------------------------------------
c  Calculate vib.-rotational constants for E_vj
c-------------------------------------------------
             call  convj(gg,300,We)
               write(31,840) De,We,w0,we0,wex,wey,wez,wet,wes,wer
	     ENDIF
C====================================================================
  534 format(//3x,'The dissociation energy    De  = ',1PE16.8,
     #//,'   Equili. intern. distance   Re  = ',1PE16.8,
     #//,'   Rydberg exponental para. alpha = ',1PE16.8,
     # /,'     alpha = We * dsqrt( amu/De ) ', 
     #//,'    Morse  exponental para. beta0 = ',1PE16.8,
     # /,'     beta0 = We * dsqrt( amu/(2*De) ) ',/
     #//,'   f_(n) are from ANAlytical derivative code. ',
C    #//,'   f_(n) are from numerical derivative code. ',
     #//9x,'  err are the errors of f_(n).',
     #///8x,'   The nth force constants are : ',
     #//5x,"  n   f_(n; Har/ao**n)     err(f) "/)
  536 format(///4x,'The nth force constants in other units are : ',
C    #//3x,' e_(n) are the force constatns from analytical ',
     #//3x,'    The force constatns are from analytical ',
     # /3x,'       derivatives of Rydberg potential. ',
     #//3x,'      1 aJ = 1 attojoule = 0.229371 Har. ',
     # /3x,'             1 ao = 0.529177249 A ',
     #//5x,"  n   f_(n; Har/ao**n)  f_(n; aJ/A**n) "/)
C    #//5x,"  n    f_(n; aJ/A**n)   e_(n; aJ/A**n) "/)
c    #//5x,"  n    f_(n; aJ/A**n)      err(f) "/)
  540 format(//3x,'The INItial GUESSES of VIB. FORCE constants : ',
     # /3x,'For  De =',f12.8,' ;   We =',1pe13.5,
     #//4x,"The potential variational parameter  b  is",/17x,1pe13.5,
     #//5x,"  n   f_(n; Har/ao**n)  f_(n; aJ/A**n) "/)
  550 format(5x,i3,2x,3(1PE16.8,x) )
  552 format(2x,i2,3x,1PE28.20)
  560 format(1PE24.16)
  730 format(///6x,'--%&*  Checking the qualities of force',
     #' constants  *&%-- ',/
     #/8x,'  As R --> inf == Rmax, physically CORRECT force ',
     #/8x,'         constants f_n should satisfy : ',
     #//8x,'V(R->inf=Rmax)_i = De = Sum_n  f_n*(Rm_i - Re)**n/n! ',
     #//7x,'    De == De_g :   Rm_g = Rmax + Rless ;  i == g ',
     # /7x,'    De == De_l :   Rm_l = Rmax - Rless ;  i == l ',
     #//10x,'For De_i(m) = De_i(Morse) :  f_n == f_n(Morse)  ',
     # /10x,'For De_i(r) = De_i(Rydbg) :  f_n == f_n(Rydbg)  ',
     #//10x,'Dconv - Tolerance used to check if V(Rmax) = De ',
     # /10x," Nryd - No of f_n's used to cal. V(Rmax) ",
     #//8x,'De',11x,'Rmax',10x,'Rless',11x,'Dconv',8x,'Nryd',
     #/1x,4(1pe14.6,x),3x,i3,
     #///'  n',5x,'De_g(m)',7x,'De_l(m)',7x,'De_g(r)',
     #7x,'De_l(r)',4x,'{De_g-De_l}_r ',/)
  732 format(36x,'V_rydberg(R) converged ',/)
  734 format( 9x,'V_morse(R) converged ',/)
  740 format(i3,x,5(1pe14.6) )
  742 format(//17x,'V(R->inf=Rmax) = De  is reached with : ',
     #//23x,'  n',5x,'De(Rydbg; a.u)',/23x,i3,3x,1pe16.8,
     #//23x,'  n',5x,'De(Morse; a.u)',/23x,i3,3x,1pe16.8,/)
  745 format(//22x,'Input We         Corrected We',
     #/22x,8(1H-),9x,12(1H-),/17x,1PE16.8,3x,1PE16.8)
  750 format(///3x,'Force constants ( h ) are from numerical ',
     #'derivative code  "dfridr" ;',
     #/3x,'Force constants ( p & f ) are from ANAlytical formulae :',
     #/3x,'         For  V_p-g(R),  2 =< n =< 11  ONLY  . ',
     #//5x,"From    V_morse(R)      V_rydberg(R)       V_p-g(R)  ",
     # /5x,"  n   p(n; Har/ao**n)  f(n; Har/ao**n)  h(n; Har/ao**n) "/)
  752 format(///3x,'Force constants ( h ) are from numerical ',
     #'derivative code  "dfridr" ;',
     #/3x,'Force constants ( p & f ) are from ANAlytical formulae :',
     #/3x,'         For  V_p-g(R),  2 =< n =< 11  ONLY  . ',/
     #/3x,'         1 aJ = 1 attojoule = 0.229371 Har; ',
     #/3x,'              1 ao = 0.529177249 A . ',
     #//5x,"From    V_morse(R)       V_rydberg(R)      V_p-g(R)  ",
     # /5x,"  n    p(n; aJ/A**n)    f(n; aJ/A**n)    h(n; aJ/A**n) "/)
  760 format(//12x,'The ERRORS for these force constants are : ',
     #//5x,"  n     e( V_morse )    e( V_rydberg )     e( V_p-g ) "/)
  790 format(///13x,'=#=  OUTPUT from "VRener.f"  =#= ',/
     #/6x,'The switches for NUMerical derivatives & E_vj : ',/
     #/2x,"[   Mwe = 1, New  We = We + bryd ;  = 0, We = We . ",/
     #/12x,"Use following V(R) to cal. n_th force constants f_n's.",
     #/8x,'|',/8x,'| = 1, use Vmorse as V(R).',
     #/8x,'|',/8x,'| = 2, use Vryd(R) = -De *[1 + a*x ]*exp(-a*x)',
     #16x,'1.)',/8x,'|',16x,
     #'a = b * dsqrt( amu/De ) ;   x = R - Re.',/8x,'|',16x,
     #'b = We + bryd ;  originally,  b = We.',
     #/8x,'|',/8x,'|',52(1H=),
     #/8x,'| The GENERALIZED Rydberg potentials are defined as : ',
     #/8x,'|',/8x,'| = 3, use',/8x,'|',6x,
     #'Uryd(R) = -Db *[exp(-a*x)/ba + (a+bryd)*x] *exp(-a*x)  2.) ',
     #/8x,'|',9x,"ba=beta;  Db=De*ba; a=We*dsqrt(amu/De) ", 
     #/8x,'|',/8x,'|',18(1H-),/8x,'| = 4, use',/8x,'|',6x,
     #'Uryd(R)= -De*[1.0 + (a-bryd)*x]*exp(-a*x)',14x,'3.) ',
     #/8x,'|',11x,'a = bryd + dsqrt(bryd*bryd + f2/De) ',
     #/8x,'|',/8x,'|',5(1H-),/8x,'| = 5, use',/8x,'|',6x,
     #'Uryd(R)= -De*(1.0 + a*x)*exp[-(a-bryd)*x]',14x,'4.) ',
     #/8x,'|',11x,'a = dsqrt(bryd*bryd + f2/De) ',
     #/8x,'|',/8x,'|',18(1H-),/8x,'| = 6, use',/8x,'|',6x,
     #'Umos(R)=De*[exp(-2*B*x) - (2+bryd*x)*exp(-B*x)]      ',
     #/8x,'|',5x,'Vcomb(R)=0.5*[ Ryd(R) + Umos(R) ]',23x,'5.) ',
     #/8x,'|',11x,'B = betam = We*dsqrt{ amu/(2*De) } ', )
  800 format(8x,'|',/8x,'|',18(1H-),/8x,'| = 7, use',/8x,'|',6x,
     #'Vmos(R)= De*[exp(-2*B*x) - 2*exp(-B*x)]      ',
     #/8x,'|',6x,'Vryd(R)=-De*[1.0 + a*x]*exp(-a*x) ',
     #/8x,'|',6x,'V(R)=Vryd(R) + b0*[ Vryd(R) - Vmos(R) ]',16x,'6.)',
     #/8x,'|',5(1H-),/8x,'| = 8, use',
     #/8x,'|',6x,'V(R)=Vmos(R) + b0*[ Vmos(R) - Vryd(R) ]',16x,'7.)',
     #/8x,'|',
     #/8x,'|',8x,'where  b0=bryd**Mryd  or  b0=(bryd/Re)**Mryd  ',
     #/8x,'|',10x,'and [...] has the same sign as the 1st term. ',
     #/8x,'|',/3x,'Ny = |',18(1H-),/8x,'| = 9, use',
     #/8x,'|',6x,'Va(R)=Vcomb(R) + b0*[ Vryd(R) - Vmos(R) ]',
     #13x,' 8.)', 
     #/8x,'|',6(1H-),/8x,'| = 10, use',
     #/8x,'|',6x,'Vb(R)=Vcomb(R) + b0*[ Vmos(R) - Vryd(R) ]',
     #13x,' 9.)',/8x,'|',/8x,'|',18(1H-),/8x,'| = 11, use',
     #/8x,'|',6x,'U(R) = [ Um(R) + Ur(R) ]/2.0 ',25x,'10.)',
     #/8x,'|  where ',
     #/8x,'|',4x,'Um(R) = De*[exp(-2*bt*x) - (2+b*x)*exp(-bt*x)] ',
     #/8x,'|',4x,'Ur(R) = -De*[1 + (af-b)*x]*exp(-af*x) ',
     #/8x,'|',4x,'bt = 0.5*(-cb + dsqrt(cb**2 + 2.0*ff2/De) ) ',
     #/8x,'|',4x,'cb = 0.5*( 1.0-sqrt(2.0) )*b ;  b == bryd',/8x,'|',
     #/8x,'|',6(1H-),/8x,'| = 12,  SAME as Eq.(10) BUT ',
     #/8x,'|',4x,'bt = 0.5*(-cb + dsqrt(cb**2 + 4.0*ff2*bk/De) )/bk',
     #/8x,'|',4x,'cb = 0.5*( 1.0 - ak)*b ;      b == bryd  ; ',
     #/8x,'|',4x,'bk = 1.0 + ak*ak/2.0 ;        af = ak*bt . ',
     #/8x,'|',4x,'ak = INPUT value closing to sqrt(2.0)=1.4142 . ',
     #/8x,'|',/8x,'|',6(1H-),/8x,'| = 13,  SAME as Eq.(10) BUT ',
     #/8x,'|',4x,'Um(R) = De*[(1-b*x)*exp(-2*bt*x) - 2*exp(-bt*x)] ',
     #/8x,'|',4x,'bt = 0.5*(-cb + dsqrt(cb**2 + 2.0*ff2/De) ) ',
     #/8x,'|',4x,'cb = ( 1.0-sqrt(2.0)/2 )*b ;  b == bryd . ',
     #/8x,'|',/8x,'|',6(1H-),
     #/8x,"| = 18, use Eq.10) & PVM procedures to calculate the ",
     #/8x,'|',7x,"INItial guesses of the FORCE constants fn's ",
     #/8x,'|',7x,"by using VIBrational constants from AM: A*X = E. ",
     #/8x,'|',/8x,'|',37(1H=),
     #/8x,'| = 26, use Vpseudo_gaussian  as V(R). ',/
     #/4x,"  Ner = No of f_n's used to cal. ro-vib constants. ",/
     #/4x,'   Ny = 4, 5, 12, 13, 18  are GOOD choices ! ',
     #/4x,'      = 4, 5    --> f_1 =/= 0.0 which is NOT physical; ',
     #/4x,'      =12,13,18 --> f_1 === 0.0 which is CORRECT ! ',/
     #/4x,"   ak - a value close to 1.4142 for Ny = 12 & used in ",
     #/4x,"          U(R) = [Um(R) + Ur(R)]/2.0 ;   af=ak*beta. ",
     #/4x," bryd - The variational constant used to adjust Vryd(R), ",
     #/4x," Mryd = M = powers in Eq 7.) for Ny = 8 . ",
     #/6x,"    bryd, Mryd are meaningless for Ny = 1 & 26 .",12x,"]",
     #//10x,'Mwe  Ner  Ny    ak         bryd         Mryd ',/
     #10x,i2,3x,i2,3x,i2,2x,f6.4,1x,f16.12,4x,i2//)
  810 format(10x,'Mwe, Ny, ak, bryd, Mryd  can be ANY value for ',
     #'mtyp = 3',/10x,23(1H*),9x,3(1H*),11x,4(1H*),//
     #//4x,'Switch for first GUESS of vibrational FORCE ',
     #"constants  f_n's : ",
     #/16x,"ksh = 0, use fn's from INPUT value above; ",
     #/16x,"    = 1, use fn's from Morse V(R);  ",
     #/16x,"    = 2, use fn's from Rydberg-like V(R);  ",
     #/16x,"    = 3, use fn's from P-G V(R).  ",/
     #/29x,"ksh =",i2,/
     #/3x,'[ hs1, hs2 --, Estimated initial stepsize ',
     #'used by the ',/20x,'numerical derivative code "dfirdr" ; ',
     #/3x,'                 hs1 for Vrydberg,   hs2 for Vp_g  .   ]',
     #//16x,'hs1 =',f8.4,4x,'hs2 =',f8.4,//)
  812 format(/2x,"Use following V(R) to cal. INItial GUESS for ",
     #"n_th FORCE constants f_n's.",/
     #/4x,'   Ny = The TYPE of Rydberg potential Vryd(R; We).',
     #/4x,'      = 4, 5, 12, 13, 18  are GOOD choices ! ',
     #/4x,'      = 4, 5    --> f_1 =/= 0.0 which is NOT physical; ',
     #/4x,'      =12,13,18 --> f_1 === 0.0 which is CORRECT ! ',/
     #/5x,"bryd = The variational constant used to adjust Vryd(R;We),",
     #//20x,"Ny             bryd ",/20x,i2,7x,f16.12,/
     #/10x," !  bryd  MAY affect the VALUES of f_n's  !  ",/) 
  820 format(//3x,'Switches to scale SOME calc. VIB-ROT constants :',
     #//3x,'[ aye, ..., are = 0, Zero calc. VIB. constants ',
     # /3x,'                     WeYe, WeZe, WeTe, WeSe, WeRe ; ',
     # /3x,'                = 1, Do NOT change sign of constants; ',
     # /3x,'                =-1, Change the sign of VIB. constants. ',
     #//3x,'  abe, aae, age, ',
     # /3x,'  ae3  ---  ae7 = 0, Zero calc. ROT. constants  ',
     # /3x,'                     B_e,Alpha_e,Gamma_e,Eta3,...,Eta7 ; ',
     # /3x,'                = 1, Do NOT change sign of constants; ',
     # /3x,'                =-1, Change the sign of ROT. constants. ',
     #//3x,'  ade, abt, ax2, ',
     # /3x,'  ax3  ---  ax7 = 0, Zero calc. ROT. constants  ',
     # /3x,'                     D_e,Beta_e,Xsi2,Xsi3, ...,Xsi7 ; ',
     # /3x,'                = 1, Do NOT change sign of constants; ',
     # /3x,'                =-1, Change the sign of ROT. constants. ]',
     #//8x,'        aye  aze  ate  ase  are   ',
     #/15x,f4.1,1x,f4.1,1x,f4.1,1x,f4.1,1x,f4.1,1x,
     #//3x,'        abe  aae  age  ae3  ae4  ae5  ae6  ae7  ',
     #/10x,f4.1,1x,f4.1,1x,f4.1,1x,f4.1,1x,f4.1,1x,f4.1,
     #1x,f4.1,1x,f4.1,
     #//3x,'        ade  abt  ax2  ax3  ax4  ax5  ax6  ax7  ',
     #/10x,f4.1,1x,f4.1,1x,f4.1,1x,f4.1,1x,f4.1,1x,f4.1,
     #1x,f4.1,1x,f4.1,///
     #/5x,"br4,br6,br8,br10,br12,br14 = Scaling switches for long ",
     #/5x,"range corce const f_n(long;R=Re) = d~n/(dr~n)[1/R**k]  ",
     #/5x,"which will be added into f_n(R=Re) = V_rydberg~n(R=Re) ",
     #/5x,"and are important for Wande-Waals molecules and some ",
     #/5x,"quasi-stable molecular electronic states . ",/
     #/5x,"f_n(long;R=Re) = (-1.0)**n *(k+n-1)!/[(k-1)!*R**(k+n)] ",
     #/8x,"brn = 1.0, Add f_n(long;R=Re); = 0.0, do NOT add. ",
     #//4x,"           br4   br6   br8  br10  br12  br14 ",
     #/14x,f4.1,2x,f4.1,2x,f4.1,2x,f4.1,2x,f4.1,2x,f4.1,/)
  840 format(/1x,'VIBrational SPECTRUM constants from ECM  ',
     #//3x," For  De = ",1pe16.8,"  a.u.",
     #//3x," We(inp) = ",1pe16.8,"  a.u.",
     # /3x,"    W0   = ",1pe16.8,"  a.u.",
     # /3x,"   We0   = ",1pe16.8,"  a.u.",
     # /3x,"   WeXe  = ",1pe16.8,"  a.u.",
     # /3x,"   WeYe  = ",1pe16.8,"  a.u.",
     # /3x,"   WeZe  = ",1pe16.8,"  a.u.",
     # /3x,"   WeTe  = ",1pe16.8,"  a.u.",
     # /3x,"   WeSe  = ",1pe16.8,"  a.u.",
     # /3x,"   WeRe  = ",1pe16.8,"  a.u." )
c----------------------------------------------------------
  900   return
      end
C===
      function dfridr(n,ns,x,h,err,We1)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
c----------------------------------------------------------
c   Returns the derivative of a function func at a point
c x by Ridders' method of polynomial extrapolation. The 
c value h is input as an estimated initial stepsize; it
c needs not be small, but rather should be an increment
c in x over which func changes substantially. An estimate
c of the error in the derivative is returned as err.
c   Parameters:  Stepsize is decreased by CON at each
c iteration. Max size of tableau is set by NTAB. Return
c when error is SAFE worse than the best so far.
c   n -- The order of derivative.
c
c             Experiments on parameters :
c          --------------------------------
c           SAFE               Results   err          h
c ---------------------------  -------  ------  -----------
c 1.001, 1.01, 1.5, 2.0, 5.0    GOOD    1*E-16  0.001 - 4.0
c
c   Results are NOT sensitive to  CON, NTAB.
c     EXTERNAL func
c-
c     PARAMETER (CON=1.4,CON2=CON*CON,BIG=1.E30,NTAB=10,SAFE=2.)
c     dimension  a(NTAB,NTAB)
c----------------------------------------------------------
      EXTERNAL fpot
      PARAMETER (NTAB=100)
      dimension  a(NTAB,NTAB)
c---------------------------------
c  ns=2 for V_rydberg(R)
c---------------------------------
      if (ns .eq. 2) then
         CON = 1.40d0
         BIG = 1.0E+30
        SAFE = 2.0
c---------------------------------
c  ns=3 for V_pseudo-gaussian(R)
c---------------------------------
      elseif (ns .eq. 3) then
         CON = 5.0d0
         BIG = 1.0E+30
        SAFE = 3.0
      endif
        CON2 = CON*CON
c----------------------------------------------------------
      if(h.eq.0.) pause 'h must be nonzero in dfridr'
        hh=h
c     a(1,1)=(func(x+hh)-func(x-hh))/(2.0*hh)
      a(1,1)=(fpot(x+hh,n,ns,We1)-fpot(x-hh,n,ns,We1))/(2.0*hh)
        err=BIG
c----------------------------------------------------------
c   Successive columns in the Neville tableau will go to
c smaller stepsizes and higher orders of extrapolation.
c----------------------------------------------------------
      do 12 i=2,NTAB
c---------------------------------
c Try new, smaller stepsize.
c---------------------------------
          hh=hh/CON
c         a(1,i)=(func(x+hh)-func(x-hh))/(2.0*hh)
          a(1,i)=(fpot(x+hh,n,ns,We1)-fpot(x-hh,n,ns,We1))/(2.0*hh)
          fac=CON2
c-------------------------------------------
c Compute extrapolations of various orders,
c requiring no new function evaluations.
c-------------------------------------------
        do 11 j=2,i
          a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1.)
          fac=CON2*fac
          errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
          if (errt.le.err) then
            err=errt
            dfridr=a(j,i)
          endif
  11    continue
          if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err)return
  12  continue
        return
      END
C
      subroutine Vpotdata
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectrb/ Be,alphae,gamae,Der,betae,bryd,wvib0(40)
      common /pgdata1/ bt,c0,c02,c03,c04,c05,c06,c07,c08,c09,c010
      common /pgdata2/ Re2,Re4,Re6,Re8,Re10,Re12,Re14,Re16,Re18,Re20
      common /vnumpot/ ff2, betap
c----------------------------------------------------------
       ff2 = amu*We*We
      betap= dsqrt( ff2/(2.0*De) )
c----------------------------------------------------------
        bt = 0.50d0*dsqrt(4.0d0 + ff2*Re*Re/De) - 1.0
        c0 = - 2.0d0*bt
c
       Re2 = Re*Re
       Re4 = Re2*Re2
       Re6 = Re4*Re2
       Re8 = Re6*Re2
      Re10 = Re8*Re2
      Re12 = Re10*Re2
      Re14 = Re12*Re2
      Re16 = Re14*Re2
      Re18 = Re16*Re2
      Re20 = Re18*Re2
c
      c02  = c0*c0
      c03  = c02*c0
      c04  = c03*c0
      c05  = c04*c0
      c06  = c05*c0
      c07  = c06*c0
      c08  = c07*c0
      c09  = c08*c0
      c010 = c09*c0
c
        return
      END
C
      function fmorse(n)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /vnumpot/ ff2, betap
c----------------------------------------------------------
c     ff2 = amu*We*We  & betap=dsqrt(ff2/(2.0*De) )
c----------------------------------------------------------
      fmorse = (-1)**n * ( 2.0**(n-1) - 1.0 ) * betap**(n-2) * ff2 
c----------------------------------------------------------
        return
      END
c
      function fryd(n,R,We1)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectrb/ Be,alphae,gamae,Der,betae,bryd,wvib0(40)
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /bdt/ ak,beta,br4,br6,br8,br10,br12,br14
c----------------------------------------------------------
c   The nth derivative of Rydberg potentials;  n >= 2. 
c----------------------------------------------------------
c        a =  We1 * dsqrt(amu/De) - bryd
         a = ( We1 + bryd ) * dsqrt(amu/De)
        ep = dexp( -a*beta*(R-Re) )
      if (Ny .eq. 2) then
        Vr = - De*beta*( 1.0/beta + a*(R-Re) )*ep
        fryd = (-1.0)**n * (a*beta)**n * ( n*De*ep + Vr )
          goto 50
c----------------------------------------------------------
c   The nth derivative of the GENERALIZED Rydberg 
c potentials.
c----------------------------------------------------------
      elseif (Ny .ge. 3 .and. Ny .le. 25) then
           a = We1 * dsqrt( amu/De )
           x = R - Re
          ep = dexp(-a*x)
          ep2= dexp(-2.0*a*x)
          a1 = (-a)**n
          a2 = (-a)**(n-1)
          Vr = - De*beta*( ep/beta - (a + bryd)*x )*ep
c-
        if (Ny .eq. 3) then
            aa = 2.0**n - 1
          fryd = a1*Vr - n*a2*De*beta*(a+bryd)*ep + aa*a1*De*ep2 
               goto 50
        endif
C-----------------------------------------------------------
C  Ny=4  for  Ur(R)=-De*[1 + (a - bryd)*x]*exp(-a*x)
C-----------------------------------------------------------
        if (Ny .eq. 4) then
             f2a =  We1 * We1*amu
               a = bryd + dsqrt(bryd*bryd + f2a/De)
             epa = dexp(-a*(R-Re))
      fryd=(-a)**n+n*(a-bryd)*(-a)**(n-1)+(-a)**n*(a-bryd)*(R-Re)   
              fryd = -fryd*De*epa
               goto 50
        endif 
C-----------------------------------------------------------
C  Ny=5  for  Ur(R)=-De*[1 + a*x] * exp[-(a-bryd)*x]
C-----------------------------------------------------------
        if (Ny .eq. 5) then
             f2a =  We1 * We1*amu
               a =  dsqrt(bryd*bryd + f2a/De)
             epa = dexp(-(a-bryd)*(R-Re))
      fryd=(-(a-bryd))**n+n*(-(a-bryd))**(n-1)*a 
      fryd= fryd + a*(-(a-bryd))**n*(R-Re)   
              fryd = -fryd*De*epa
               goto 50
        endif 
C-----------------------------------------------------------
            Ur = - De*( 1.0 + (a - bryd)*x )*ep
          temp = a1*Ur - n*a2*(a - bryd)*De*ep
           b = We1 * dsqrt( amu/(2.0*De) )
          eb = dexp(-b*x)
          eb2= dexp(-2.0*b*x)
          b1 = (-b)**n
          b2 = (-b)**(n-1)
          b3 = (-2.0*b)**n
          tem0 = b3*eb2 - n*bryd*b2*eb - b1*(2.0 + bryd*x)*eb 
          vcomb = (De*tem0 + temp)/2.0
        if (Ny .eq. 6) then
          fryd = vcomb
               goto 50
        endif
c-
          vrn = a1*( Vr + n*De*ep)
          vmn = 2.0*De*( 2**(n-1) *eb2 - eb ) * (-b)**n
c
          bm = (bryd/Re)**Mryd          
C         bm = bryd**Mryd          
             b0 = bm
C            b0 = bryd
          if (Ny .eq. 7) then
	        difv = vrn - vmn
            if ( vrn .gt. 0.0 ) then
	          fryd = vrn + b0*difv
		  if ( vrn .lt. vmn) fryd = vrn - b0*difv
            else
	          fryd = vrn - b0*difv
		  if ( vrn .lt. vmn) fryd = vrn + b0*difv
            endif
               goto 50
          endif
c
          if (Ny .eq. 8) then
	        difv = vmn - vrn
            if ( vmn .gt. 0.0 ) then
	          fryd = vmn + b0*difv
		  if ( vmn .lt. vrn) fryd = vmn - b0*difv
            else
	          fryd = vmn - b0*difv
		  if ( vmn .lt. vrn) fryd = vmn + b0*difv
            endif
               goto 50
          endif
c
          if (Ny .eq. 9) then
	         difv = vrn - vmn
            if ( vcomb .gt. 0.0 ) then
	          fryd = vcomb + b0*difv
		  if ( vrn .lt. vmn) fryd = vcomb - b0*difv
            else
	          fryd = vcomb - b0*difv
		  if ( vrn .lt. vmn) fryd = vcomb + b0*difv
            endif
               goto 50
          endif
c
          if (Ny .eq. 10) then
	         difv = vmn - vrn
            if ( vcomb .gt. 0.0 ) then
	          fryd = vcomb + b0*difv
		  if ( vmn .lt. vrn) fryd = vcomb - b0*difv
            else
	          fryd = vcomb - b0*difv
		  if ( vmn .lt. vrn) fryd = vcomb + b0*difv
            endif
               goto 50
          endif
C-----------------------------------------------------------
C       Ny=11; 13; 18  for  U(R) = [ Um(R) + Ur(R) ]/2.0  
C  where 
C     Ur(R) = -De*[1 + (af-b)*x]*exp(-af*x)
C  bt = 0.5*( -ccb*b + dsqrt((ccb*b)**2 + 2.0*ff2/De))
C-
C        For Ny = 11 :
C     Um(R) = De*[exp(-2*bt*x) - (2+b*x)*exp(-bt*x)]
C       ccb = 0.5*( 1.0-sqrt(2.0) );  b == bryd.
C-
C        For Ny = 13; 18 :
C     Um(R) = De*[(1-b*x)*exp(-2*bt*x) - 2*exp(-bt*x)]
C       ccb = 1.0-sqrt(2.0)/2 ;  b == bryd.
C-----------------------------------------------------------
      if (Ny .eq. 11 .or. Ny .eq. 13 .or. Ny .eq. 18) then
        ff2 = We1*We1*amu
        sqr2= dsqrt(dfloat(2.0))
          if (Ny .eq. 11) ccb = 0.5*(1.0-sqr2)
          if (Ny .eq. 13 .or. Ny .eq. 18) ccb = 1.0 - sqr2/2.0
         bt = 0.5*( -ccb*bryd + dsqrt((ccb*bryd)**2 + 2.0*ff2/De))
         af = sqr2*bt
        ep2 = dexp( -2.0*bt*(R-Re))
        ep1 = dexp(     -bt*(R-Re))
        epa = dexp(     -af*(R-Re))
        Uryd= -De*( 1.0 + (af - bryd)*(R - Re) )*epa
C
        if (Ny .eq. 11) then
          fryd = (-2.0*bt)**n *ep2 - n*bryd*(-bt)**(n-1) *ep1
          fryd = fryd - (-bt)**n *(2.0+bryd*(R-Re))*ep1
        elseif (Ny .eq. 13 .or. Ny .eq. 18) then
          fryd = ( 1.0 - bryd*(R - Re) )*(-2.0*bt)**n *ep2 
          fryd = fryd - n*bryd*(-2.0*bt)**(n-1) *ep2
          fryd = fryd - 2.0*(-bt)**n * ep1
        endif
          fryd = De*fryd
C End Morse Function  and add Rydberg.
          fryd = fryd + (-af)**n *Uryd
          fryd = fryd - n*(-af)**(n-1) *(af - bryd)*De*epa
          fryd = 0.5*fryd
            goto 50
      endif
C-----------------------------------------------------------
C         Ny=12 for  U(R) = [ Um(R) + Ur(R) ]/2.0  
C  where 
C     Um(R) = De*[exp(-2*bt*x) - (2+b*x)*exp(-bt*x)]
C     Ur(R) = -De*[1 + (af-b)*x]*exp(-af*x)
C        bt =0.5*(-ccb + dsqrt(ccb**2 + 4.0*ff2*bk/De) )/bk
C       ccb = (1.0 - ak)*b;  b == bryd ;  ak = INPUT value ;
C        bk = 1.0 + ak*ak/2.0 ;      af = ak*bt .
C-----------------------------------------------------------
      if (Ny .eq. 12 ) then
        ff2 = We1*We1*amu
        ccb = (1.0 - ak)*bryd
         bk = 1.0 + ak*ak/2.0
         bt =0.5*(-ccb + dsqrt(ccb**2 + 4.0*ff2*bk/De) )/bk
         af = ak*bt
        ep2 = dexp( -2.0*bt*(R-Re))
        ep1 = dexp(     -bt*(R-Re))
        epa = dexp(     -af*(R-Re))
        Uryd= -De*( 1.0 + (af - bryd)*(R - Re) )*epa
C
       fryd = (-2.0*bt)**n *ep2 - n*bryd*(-bt)**(n-1) *ep1
       fryd = fryd - (-bt)**n *(2.0+bryd*(R-Re))*ep1
       fryd = De*fryd
C End Morse Function  and begin Rydberg.
        fryd = fryd + (-af)**n *Uryd
        fryd = fryd - n*(-af)**(n-1) *(af - bryd)*De*epa
        fryd = 0.5*fryd
          goto 50
      endif
C-----------------------------------------------------------
C  
          if (Ny.eq.14 ) then
C     fryd = Rdv(n,R,We1)   
      fryd = RMdv(n,R)   
          do i=0,n
      fryd = fryd + facx(n)/(facx(n-i)*facx(i))
c     fryd = fryd - facx(n)/(facx(n-i)*facx(i))
     #           *Gdv((n-i),R)*(Rdv(i,R,We1)-RMdv(i,R))  
          enddo     
              goto 50
          endif 
c---
          if ( Ny .eq. 15 ) then
         a = We1 * dsqrt( amu/De )
        ah = ( We1 + bryd ) * dsqrt( amu/De )
       ep  = dexp ( -a * (R - Re) ) 
      fryd0 = ( -a )**n
      fryd  = 0.0            
      do m = 1 ,int(Mryd)
         if ( n .ge. 2*m -1 )then
      fryd = fryd + Cnj(n,2*m-1)*(-a)**(n-2*m+1) * ah**(2*m-1)   
c     fryd = fryd + Cnj(n,2*m-1)*(-a)**(n-2*m+1)*facx(2*m-1)
c     fryd = fryd + Cnj(n,j)*(-a)**(n-j)*(R-Re)**(2*m-1-j)*Pnj(2*m-1,j)
         endif
      enddo
      fryd = (fryd + fryd0) * ( - De )      
            goto 50
          endif
C----
          if ( Ny .eq. 16 ) then
         a = ( We1 + bryd ) * dsqrt( amu/De )
        ah = We1 * dsqrt( amu/De )
       ep  = dexp ( -a * (R - Re) ) 
        fryd = 0.0
        fryd0= ( -a )**n
      do m = 1 , int(Mryd)
        if ( n .ge. 2*m -1 )then
          fryd = fryd + Cnj(n,2*m-1)*(-a)**(n-2*m+1) * ah**(2*m-1)
c-
c       kk = 2*m - 1
c     do j = 0 , 2*m - 1    
c     do j = 0 , kk         
c     fryd = fryd + Cnj(n,j)*(-a)**(n-j)*(R-Re)**(2*m-1-j)*Pnj(2*m-1,j)
c        if (n.ge.j) then
c     fryd = fryd + facx(n)/(facx(j)*facx(n-j))*(-a)**(n-j)
c    # *(R-Re)**(kk-j) *facx(n)/facx(n-j) * ah**(2*m-1)/facx(2*m-1)
c        else
c-
        endif
      enddo
c-
c       fryd = fryd * ah**(2*m-1)/facx(2*m-1)
c     enddo
c-
      fryd = fryd + fryd0               
      fryd = fryd * ep * ( - De )      
            goto 50
          endif
C----
      if ( Ny .eq. 17 ) then
c        a = We1 * dsqrt( amu/De )
         ff2b= We1*We1*amu 
         a =  bryd + dsqrt( bryd*bryd + ff2b/De)
C1       a = -bryd + dsqrt( bryd*bryd + ff2b/De)
c       ah = ( We1 + bryd ) * dsqrt( amu/De )
C1      ah =  a + bryd
c-
        ah =  a - bryd
       ep  = dexp ( -a * (R - Re) ) 
      fryd0 = ( -a )**n
      fryd  = 0.0            
c-
      do m = 1 ,int(Mryd)
         if ( n .ge. 2*m -1 )then
      fryd = fryd + Cnj(n,2*m-1)*(-a)**(n-2*m+1) * ah**(2*m-1)   
c     fryd = fryd + Cnj(n,2*m-1)*(-a)**(n-2*m+1)*facx(2*m-1)
c     fryd = fryd + Cnj(n,j)*(-a)**(n-j)*(R-Re)**(2*m-1-j)*Pnj(2*m-1,j)
         endif
      enddo
        fryd = (fryd + fryd0) * ( - De )      
            goto 50
      endif
C
C-----------------------------------------------------------
C  Ny=19 for  Ur(R)=
C-----------------------------------------------------------
      if ( Ny .eq. 19 ) then
        ff2 = We1*We1*amu
         bt = (-bryd + dsqrt(bryd*bryd + 2.0*ff2/De))/2.0
         fryd = 0.0
        ep2 = dexp(-2.0*bt*(R-Re))
        ep1 = dexp( -bt*(R-Re)   )
        fryd = (-2.0*bt)**n*ep2 - n*bryd*(-bt)**(n-1.0)*ep1
        fryd = fryd - (-bt)**n*(2.0+bryd*(R-Re))*ep1
        fryd = De*fryd
          goto 50
      endif
C-----------------------------------------------------------
        if (Ny .eq. 20 ) then
          ff2 = We1*We1*amu
          sqr2= dsqrt(Dfloat(2.0))
          ccb = 0.5*(1.0+sqr2)
         bt = 0.5*( -ccb*bryd + dsqrt((ccb*bryd)**2 + 2.0*ff2/De))
           af = sqr2*bt
          ep2 = dexp( -2.0*bt*(R-Re))
          ep1 = dexp(     -bt*(R-Re))
          epa = dexp(     -af*(R-Re))
          Uryd= -De*(1 + (af + bryd)*(R - Re))*epa
C End Morse Function  and begin Rydberg.
          fryd = (-2.0*bt)**n*ep2 - n*bryd*(-bt)**(n-1)*ep1
          fryd = fryd - (-bt)**n*(2.0+bryd*(R-Re))*ep1
          fryd = De*fryd
C
          fryd = fryd + (-af)**n*Uryd
          fryd = fryd - n*(-af)**(n-1)*(af + bryd)*De*epa
          fryd = 0.5*fryd
            goto 50
        endif
C===
      endif
c--------------------------------------------------------
c   Add  d~n/(dR~n)[1/R**n] terms;  n = 4, 6, ..., 14.
c These terms are important for Wande-Waals molecules
c and quasi-stable molecules.
c f_n(long;R=Re) = (-1.0)**n *(k+n-1)!/[(k-1)!*R**(k+n)]
c--------------------------------------------------------
 50      R4n = 0.0
         R6n = 0.0
         R8n = 0.0
        R10n = 0.0
        R12n = 0.0
        R14n = 0.0
        temp = 0.0
        pref = (-1.0)**n
          e4 = 6.0
          e6 = 120.0
          e8 = 5040.0
         e10 = 362880.0
         e12 = 39916800.0
         e14 = 6227020800.0
      if ( br4 .gt. 0.0)  R4n = pref*facx(3+n)/( e4*R**(4+n) )    
      if ( br6 .gt. 0.0)  R6n = pref*facx(5+n)/( e6*R**(6+n) )    
      if ( br8 .gt. 0.0)  R8n = pref*facx(7+n)/( e8*R**(8+n) )    
      if (br10 .gt. 0.0) R10n = pref*facx(9+n)/( e10*R**(10+n) )    
      if (br12 .gt. 0.0) R12n = pref*facx(11+n)/( e12*R**(12+n) )    
      if (br14 .gt. 0.0) R14n = pref*facx(13+n)/( e14*R**(14+n) )    
          temp = br4*R4n + br6*R6n + br8*R8n + br10*R10n
          temp = temp + br12*R12n + br14*R14n
          fryd = fryd + bryd*temp
c--------------------------------------------------------
        return
      END
c
      function Gdv(n,R)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectrb/ Be,alphae,gamae,Der,betae,bryd,wvib0(40)
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /bdt/ ak,beta,br4,br6,br8,br10,br12,br14
C
       Gdv = bryd*(-bryd**Mryd/Re)**n*Dexp(-bryd**Mryd*(R-Re)/Re)
c      Gdv = bryd*(-bryd**Mryd)**n*Dexp(-bryd**Mryd*(R-Re))
      end  
C==
      function Rdv(n,R,We1)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectrb/ Be,alphae,gamae,Der,betae,bryd,wvib0(40)
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /bdt/ ak,beta,br4,br6,br8,br10,br12,br14
c----------------------------------------------------------
c   The nth derivative of Rydberg potentials;  n >= 2. 
c----------------------------------------------------------
c        a = ( We1 + bryd ) * dsqrt(amu/De)
         a =  We1 * dsqrt(amu/De)
        ep = dexp( -a*(R-Re) )
        Vr = - De*( 1.0 + a*(R-Re) )*ep
        Rdv = (-1.0)**n * a**n * ( n*De*ep + Vr )
c----------------------------------------------------------
         return
      end      


      function  RMdv(n,R)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /vnumpot/ ff2, betap
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
c----------------------------------------------------------
      if(n.eq.0) then
       RMdv=De*(dexp(-betap*(R-Re))-2.0)*dexp(-betap*(R-Re))
      else
      RMdv= (-1)**n * ( 2.0**(n-1) - 1.0 ) * betap**(n-2) * ff2 
      endif
c----------------------------------------------------------
         return
      END
c===
      function fpot(R,n,ns,We1)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectrb/ Be,alphae,gamae,Der,betae,bryd,wvib0(40)
      common /pgdata1/ bt,c0,c02,c03,c04,c05,c06,c07,c08,c09,c010
      common /pgdata2/ Re2,Re4,Re6,Re8,Re10,Re12,Re14,Re16,Re18,Re20
      common /vnumpot/ ff2, betap
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /bdt/ ak,beta,br4,br6,br8,br10,br12,br14
c----------------------------------------------------------
      if (ns .eq. 2) goto  100
      if (ns .eq. 3) goto  200
c----------------------------------------------------------
c   The nth derivative of Rydberg potentials;  n >= 2.
c 	  f(2) = amu*We1*We1 ;  
c       d/(dR){ f[ V_rydberg(R) ] } = {fpot}'
c
c100    a = dsqrt(ff2/De)
c100    a = We1 * dsqrt(amu/De)
c----------------------------------------------------------
 100  if (Ny .eq. 2) then
        a = ( We1 + bryd ) * dsqrt(amu/De)
      endif
        m = n - 1
      ep = dexp( -a*beta*(R-Re) )
      Vr = - De*beta*( 1.0/beta + a*(R-Re) )*ep
      fpot = (-a*beta)**m * ( m*De*ep + Vr )
        goto  800
c----------------------------------------------------------
c   The nth derivative of V(Pseudo-Gaussian; R) ; 
c        2 =< n =< 11  ONLY .
c 	V(n)[p-g; R] = d/(dR){ V(n-1)[p-g; R] } = {fpot}'
c----------------------------------------------------------
 200  b2 = 1.0d0 - R*R/(Re*Re)
      g1 = De*dexp(bt*b2)
      g2 = bt*g1
c
      R2 = R*R
      R3 = R2*R
      R4 = R3*R
      R5 = R4*R
      R6 = R5*R
      R7 = R6*R
      R8 = R7*R
      R9 = R8*R
      R10= R9*R
      R11= R10*R
      R12= R11*R
c
        h1 = 0.0
        h2 = 0.0
      if (n .eq. 2) then
        h1 = c0*Re2/R3 - c0*R/Re2
        h2 = c0/R - c0*R/Re2
          goto  300     
      elseif (n .eq. 3) then
        h1 = c02/R2 - 3.0*c0*Re2/R4 - c0/Re2 - c02*R2/Re4
        h2 = -c0/R2 - c0/Re2 + c02/Re2 - c02*R2/Re4
          goto  300     
      elseif (n .eq. 4) then
        h1 = -3.0*c02*R/Re4 - c03*R3/Re6 + c03/(R*Re2)
        h1 = h1 - 5.0*c02/R3 + 12.0*c0*Re2/R5
        h2 = -3.0*c02*R/Re4 + c03*R/Re4 - c03*R3/Re6
        h2 = h2 - c02/(R*Re2) + 2.0*c0/R3
          goto  300     
      elseif (n .eq. 5) then
        h1 = -3.0*c02/Re4 + c04/Re4 - 6.0*c03/(R2*Re2) -27.0*c02/R4
        h1 = h1 - 60.0*c0*Re2/R6 - 6.0*c03*R2/Re6 - c04*R4/Re8
        h2 = -3.0*c02/Re4 + 3.0*c02/(R2*Re2) - 6.0*c0/R4 
        h2 = h2 - 6.0*c03*R2/Re6 + c04*R2/Re6 - c04*R4/Re8
          goto  300     
      elseif (n .eq. 6) then
        h1 = 39.0*c03/(R3*Re2) - 168.0*c02/R5 + 360.0*c0*Re2/R7
        h1 = h1 - 15.0*c03*R/Re6 - 10.0*c04*R3/Re8 + c05*R/Re6
        h1 = h1 - 6.0*c04/(R*Re4) - c05*R5/Re10
        h2 = - 12.0*c02/(R3*Re2) + 24.0*c0/R5 - 15.0*c03*R/Re6
        h2 = h2 + 2.0*c04*R/Re6 - 10.0*c04*R3/Re8 + 3.0*c03/(R*Re4)
        h2 = h2 + c05*R3/Re8 - c05*R5/Re10
          goto  300     
      elseif (n .eq. 7) then
        h1 = -285.0*c03/(R4*Re2) + 1200.0*c02/R6 - 2520.0*c0*Re2/R8 
        h1 = h1 - 15.0*c03/Re6 - 45.0*c04*R2/Re8 - 5.0*c05/Re6
        h1 = h1 + 45.0*c04/(R2*Re4) - 15.0*c05*R4/Re10 
        h1 = h1 + c06*R2/Re8 - c06*R6/Re12
        h2 =  60.0*c02/(R4*Re2) - 120.0*c0/R6 - 15.0*c03/Re6 
        h2 = h2 + 5.0*c04/Re6 - 45.0*c04*R2/Re8 - 15.0*c03/(R2*Re4)
        h2 = h2 + 5.0*c05*R2/Re8 - 15.0*c05*R4/Re10 + c06*R4/Re10 
        h2 = h2 - c06*R6/Re12 
          goto  300     
      elseif (n .eq. 8) then
        h1 = 2340.0*c03/(R5*Re2) - 9720.0*c02/R7 + 20160.0*c0*Re2/R9
        h1 = h1 - 105.0*c04*R/Re8 - 375.0*c04/(R3*Re4) 
        h1 = h1 - 105.0*c05*R3/Re10 - 21.0*c06*R5/Re12
        h1 = h1 - 3.0*c06*R/Re8 + 45.0*c05/(R*Re6) + c07*R3/Re10
        h1 = h1 - c07*R7/Re14
        h2 = -360.0*c02/(R5*Re2) + 720.0*c0/R7 - 105.0*c04*R/Re8
        h2 = h2 + 90.0*c03/(R3*Re4) + 15.0*c05*R/Re8 
        h2 = h2 - 105.0*c05*R3/Re10 + 9.0*c06*R3/Re10 
        h2 = h2 - 21.0*c06*R5/Re12 - 15.00*c04/(R*Re6) + c07*R5/Re12 
        h2 = h2 - c07*R7/Re14
          goto  300     
      elseif (n .eq. 9) then
        h1 = -21420.*c03/(R6*Re2) + 88200.0*c02/R8 
        h1 = h1 - 181440.0*c0*Re2/R10 - 105.0*c04/Re8 
        h1 = h1 + 3465.*c04/(R4*Re4) - 420.*c05*R2/Re10 
        h1 = h1 + 42.0*c06/Re8 - 210.0*c06*R4/Re12 
        h1 = h1 - 420.0*c05/(R2*Re6) - 28.0*c07*R6/Re14 
        h1 = h1 + c08*R4/Re12 - c08*R8/Re16
        h2 = 2520.0*c02/(R6*Re2) - 5040.0*c0/R8 - 105.0*c04/Re8
        h2 = h2 - 630.0*c03/(R4*Re4) - 420.0*c05*R2/Re10 
        h2 = h2 + 42.*c06*R2/Re10 - 210.0*c06*R4/Re12 
        h2 = h2 + 105.0*c04/(R2*Re6) + 14.0*c07*R4/Re12 
        h2 = h2 - 28.0*c07*R6/Re14 + c08*R6/Re14 - c08*R8/Re16
          goto  300     
      elseif (n .eq. 10) then
        h1 = 216720.0*c03/(R7*Re2) - 887040.0*c02/R9 
        h1 = h1 + 1814400*c0*Re2/R11
        h1 = h1 - 35280.0*c04/(R5*Re4) - 945.0*c05*R/Re10 
        h1 = h1 - 1260.0*c06*R3/Re12 + 4305.0*c05/(R3*Re6) 
        h1 = h1 - 378.0*c07*R5/Re14 + 4.0*c08*R3/Re12
        h1 = h1 - 36.0*c08*R7/Re16 - 420.0*c06/(R*Re8) 
        h1 = h1  + 42.0*c07*R/Re10 + c09*R5/Re14 - c09*R9/Re18
        h2 = -20160.0*c02/(R7*Re2) + 40320.*c0/R9 
        h2 = h2 + 5040.0*c03/(R5*Re4) - 945.0*c05*R/Re10
        h2 = h2 + 84.0*c06*R/Re10 - 1260.0*c06*R3/Re12
        h2 = h2 - 840.0*c04/(R3*Re6) + 98.0*c07*R3/Re12 
        h2 = h2 - 378.0*c07*R5/Re14 + 20.0*c08*R5/Re14
        h2 = h2 - 36.0*c08*R7/Re16 + 105.0*c05/(R*Re8)
        h2 = h2 + c09*R7/Re16 - c09*R9/Re18
          goto  300     
      elseif (n .eq. 11) then
        h1 = -2404080.0*c03/(R8*Re2) + 9797760.0*c02/R10 
        h1 = h1 - 19958400.0*c0*Re2/R12 + 393120.0*c04/(R6*Re4)
        h1 = h1 - 945.0*c05/Re10 - 4725.0*c06*R2/Re12 
        h1 = h1 - 48195.0*c05/(R4*Re6) - 3150.0*c07*R4/Re14
        h1 = h1 + 54.0*c08*R2/Re12 - 630.0*c08*R6/Re16 
        h1 = h1 - 378.0*c07/Re10 + 4725.0*c06/(R2*Re8)
        h1 = h1 + 9.0*c09*R4/Re14 - 45.0*c09*R8/Re18
        h1 = h1 + c010*R6/Re16 - c010*R10/Re20
        h2 = 181440.0*c02/(R8*Re2) - 362880.0*c0/R10 
        h2 = h2 - 45360.0*c03/(R6*Re4) - 945.0*c05/Re10
        h2 = h2 + 189.0*c06/Re10 - 4725.0*c06*R2/Re12
        h2 = h2 + 7560.0*c04/(R4*Re6) + 378.0*c07*R2/Re12
        h2 = h2 - 3150.0*c07*R4/Re14 + 198.0*c08*R4/Re14 
        h2 = h2 - 630.0*c08*R6/Re16 - 945.0*c05/(R2*Re8) 
        h2 = h2 + 27.0*c09*R6/Re16 - 45.0*c09*R8/Re18
        h2 = h2 + c010*R8/Re18 - c010*R10/Re20
      endif
c---
 300  fpot = g1*h1 + g2*h2
c---
c       goto  800
c----------------------------------------------------------
 800    return
      END
C
C===
      subroutine  calEvj(nm,nv,nj,lda,Lj,Ev,Evj)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectrb/ Be,alphae,gamae,Der,betae,bryd,wvib0(40)
      common /spectr0/ w0,we0,wei,wex,wey,wez,wet,wes,wer,v0,r0
      common /spectr1/ Bee,ale,gae,eta3,eta4,eta5,eta6,eta7
      common /spectr2/ Dee,bete,xsi2,xsi3,xsi4,xsi5,xsi6,xsi7
      common /spectr4/ w08,w09,w10,w11,w12,w13,w14,w15,w16
      common /Eswitch/ aye,aze,ate,ase,are
      common /Eswitc1/ abe,aae,age,ae3,ae4,ae5,ae6,ae7
      common /Eswitc2/ ade,abt,ax2,ax3,ax4,ax5,ax6,ax7
      common /Eswitc4/ cxe,cye,cze,cte,cse,cre
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /LA1/ Deinp,Ner,mset,mcst,mtd,mdp,mv,nvs(20,80)
      complex*16  DisE3ca, DisE3c, DisE3ia, DisE3i 
      dimension  Em(200),Eu(200),Ew(200),Ex(200),Eh(200)
      dimension  e1(200),e2(200),e3(200),Ew0(200),Eu0(200)
      dimension  Ev2(200),Ev3(200),Ev4(200),Ev5(200),Ev6(200)
      dimension  Ewj(200,200),Euj(200,200),Ej(200),Er(200)
      dimension  aj1(200,200),aj2(200,200)
      dimension  Evj(Lj,Lj), Ev(Lj), Eri(200), nj(Lj)
c----------------------------------------------------------
c  nv = the number of vibrational states used.
c         (nv-1) - The HIGHEST vibrational state.
c  nj - Number of rotational states in each v state.
c==========================================================
       rydev = 13.60569809d0
        auev = 27.21139618d0
        aucm = 219474.6306d0
      rinert = amu*Re*Re
      Brigid = 1.0/(2.0*rinert)
c----------------------------------------------------------
c Goto 40 to calc. DISSOCIATION energy De using INPUT Ev's
c----------------------------------------------------------
       if (nv .eq. 1 .and. nj(1) .eq. 1) goto 20
c----------------------------------------------------------
c  Write out vib-rot constants
c----------------------------------------------------------
          ii0 = 2
           ij = 6
        if (nm .gt. 0) then
	    ii0 = 1
           ij = 38
        endif
      do i=1,ii0
        if (nm .gt. 0)  goto 10
          write(ij,400)
        if (WeXe .ne. 0.) px = 100.0*(abs(wex) - abs(WeXe))/wex
        if (WeYe .ne. 0.) py = 100.0*(abs(wey) - abs(WeYe))/wey
        if (WeZe .ne. 0.) pz = 100.0*(abs(wez) - abs(WeZe))/wez
        if (WeTe .ne. 0.) pt = 100.0*(abs(wet) - abs(WeTe))/wet
        if (WeSe .ne. 0.) ps = 100.0*(abs(wes) - abs(WeSe))/wes
        if (WeRe .ne. 0.) pr = 100.0*(abs(wer) - abs(WeRe))/wer
        write(ij,410) w0,we0,We,Wee,wex,WeXe,px,wey,WeYe,py,wez,WeZe,
     #   pz,wet,WeTe,pt,wes,WeSe,ps,wer,WeRe,pr
        write(ij,420)
        write(ij,410) w0*aucm,we0*aucm,We*aucm,Wee*aucm,
     #   wex*aucm,WeXe*aucm,px,wey*aucm,WeYe*aucm,py,
     #   wez*aucm,WeZe*aucm,pz,wet*aucm,WeTe*aucm,pt,
     #   wes*aucm,WeSe*aucm,ps,wer*aucm,WeRe*aucm,pr
c===
  10    if (nm .eq. 0) write(ij,430) 
        if (nm .gt. 0) write(ij,435) 
c
          if (Be .ne. 0.0)     rb = 100.0*(abs(Bee)- abs(Be))/Bee
          if (alphae .ne. 0.0) ra = 100.0*(abs(ale)- abs(alphae))/ale
          if (gamae .ne. 0.0)  rg = 100.0*(abs(gae)- abs(gamae))/gae
          if (Der .ne. 0.0)    sd = 100.0*(abs(Dee)- abs(Der))/Dee
          if (betae .ne. 0.0)  sb = 100.0*(abs(bete)-abs(betae))/bete
        write(ij,440) Bee,Be,rb,ale,alphae,ra,gae,gamae,rg,eta3,eta4,
     #   eta5,eta6,eta7,     Dee,Der,sd,bete,betae,sb,xsi2,xsi3,xsi4,
     #   xsi5,xsi6,xsi7
        write(ij,450) 
        write(ij,440) Bee*aucm,Be*aucm,rb,ale*aucm,alphae*aucm,
     #   ra,gae*aucm,gamae*aucm,rg,eta3*aucm,eta4*aucm,
     #   eta5*aucm,eta6*aucm,eta7*aucm,   
     #   Dee*aucm,Der*aucm,sd,bete*aucm,betae*aucm,sb,xsi2*aucm,
     #   xsi3*aucm,xsi4*aucm,xsi5*aucm,xsi6*aucm,xsi7*aucm  
           if (nm .eq. 0 .and. i .eq. 1) ij = 35
      enddo
c----------------------------------------------------------
           if (nm .gt. 0) goto 900
c----------------------------------------------------------
c  Check if the vibrational constants reproduce De
c       c : calculated data;  i : input data.
c----------------------------------------------------------
c--- From quadratic form :
        DisE2c = We*We/(4.0*wex)
        DisE2i = Wee*Wee/(4.0*WeXe)
c--- From CUBIC form :
          cubcal = wex**2 - 3.0*We*wey
        IF (cubcal .ge. 0.0) THEN
          DisE3ca= (1.0, 0.0)*2.0*( dsqrt(cubcal) )**3
        ELSE
c--- 
c--- As cubcal < 0.0, dsqrt(cubcal) is complex !
c--- ( dsqrt(cubcal) )**3 = ( i * dsqrt(- cubcal) )**3
c--- ( i * dsqrt(cubcal) )**3 = -i * ( dsqrt(- cubcal) )**3
c--- 
          DisE3ca= -(0.0, 1.0)*2.0*( dsqrt(- cubcal) )**3
        ENDIF
          DisE3cb= wex*(2.0*wex**2 - 9.0*We*wey)
          DisE3c = DisE3ca - DisE3cb*(1.0, 0.0)
          DisE3c = DisE3c/(27.0*wey*wey)
c
          cubinp = WeXe**2 - 3.0*Wee*WeYe
        IF (cubinp .ge. 0.0) THEN
          DisE3ia= (1.0, 0.0)*2.0*( dsqrt(cubinp) )**3
        ELSE
          DisE3ia= -(0.0, 1.0)*2.0*( dsqrt(- cubinp) )**3
        ENDIF
          DisE3ib= WeXe*(2.0*WeXe**2 - 9.0*Wee*WeYe)
          DisE3i = DisE3ia - DisE3ib*(1.0, 0.0)
        if (abs(WeYe) .gt. 1.0E-50) then
          DisE3i = DisE3i/(27.0*WeYe*WeYe)
        else
          DisE3i = (0.0, 0.0)
        endif
C
     
        IF( mtyp.eq. 2 ) THEN
C   Hou added this switch on 4th,July,2000 
c
c--- From quadrupole form :
c
          if ( abs(wez) .gt. 0.0 .or. abs(WeZe) .gt. 0.0 ) then
            call EvmaxDe4th(delc,deli,y4c,y4i,DisE4c,DisE4i)
          else
            DisE4i = 0.0
          endif 
c
c--- From QUINTIC form :
c
          if ( abs(wet) .gt. 0.0 .or. abs(WeTe) .gt. 0.0 ) then
            call EvmaxDe5th(y5c,y5i,DisE5c,DisE5i)
          else
            DisE5i = 0.0
          endif 
C      
        ENDIF
c
c=== Calculate percentage ERRORs :
c
        IF ( Deinp .gt. 0.0 ) THEN
          E2cp = 100.0*abs( DisE2c - De )/De
          E2ip = 100.0*abs( DisE2i - De )/De
c---
          E3cp = 100.0*abs( real(DisE3c - De) )/De
          E3ip = 100.0*abs( real(DisE3i - De) )/De
c---
          E4cp = 100.0*abs( DisE4c - De )/De
          E4ip = 100.0*abs( DisE4i - De )/De
c---
          E5cp = 100.0*abs( DisE5c - De )/De
          E5ip = 100.0*abs( DisE5i - De )/De
        ENDIF
c---
        if ( WeYe .eq. 0.0 ) E3ip = 0.0
        if ( WeZe .eq. 0.0 ) E4ip = 0.0
        if ( WeTe .eq. 0.0 ) E5ip = 0.0
c---
          if (mtyp .eq. 1 .and. abs(WeYe) .eq. 0.0)  then
	      DisE3c = 0.0
	        E3cp = 0.0
          endif
          if (mtyp .eq. 1 .and. abs(WeZe) .eq. 0.0)  then
	      DisE4c = 0.0
	        E4cp = 0.0
          endif
            if (DisE4c .eq. 0.0) E4cp = 0.0
            if (DisE4i .eq. 0.0) E4ip = 0.0
          if (mtyp .eq. 1 .and. abs(WeTe) .eq. 0.0)  then
	      DisE5c = 0.0
	        E5cp = 0.0
          endif
c----------------------------------------------------------
c   Print and compare the calculated De to see the quality
c of the constants (We, wex, wey, ...)
c----------------------------------------------------------
        write(6,460)   De, 2, DisE2c, E2cp, DisE2i, E2ip, 3,
     # real(DisE3c), E3cp, real(DisE3i), E3ip, 4, DisE4c,
     # E4cp, DisE4i, E4ip, 5, DisE5c, E5cp, DisE5i, E5ip
        write(35,460)  De, 2, DisE2c, E2cp, DisE2i, E2ip, 3,
     # real(DisE3c), E3cp, real(DisE3i), E3ip, 4, DisE4c,
     # E4cp, DisE4i, E4ip, 5, DisE5c, E5cp, DisE5i, E5ip
c---
          IF (cubcal .lt. 0.0) THEN
            write(6,470) wex**2, -3.0*We*wey, 
     # DisE3cb, DisE3ca, DisE3c
            write(35,470) wex**2, -3.0*We*wey, 
     # DisE3cb, DisE3ca, DisE3c
          ENDIF
          IF (cubinp .lt. 0.0) THEN
            write(6,475) WeXe**2, -3.0*Wee*WeYe, 
     # DisE3ib, DisE3ia, DisE3i
            write(35,475) WeXe**2, -3.0*Wee*WeYe, 
     # DisE3ib, DisE3ia, DisE3i
          ENDIF
c---
          IF (delc .lt. 0.0) THEN
            write(6,480)  delc, y4c, DisE4c
            write(35,480) delc, y4c, DisE4c
          ENDIF
          IF (deli .lt. 0.0) THEN
            write(6,485)  deli, y4i, DisE4i
            write(35,485) deli, y4i, DisE4i
          ENDIF
c----------------------------------------------------------
c   Find the MINIMUM POSITIVE root of energy derivative
c using numerical method.
c----------------------------------------------------------
             nvp = nv
         call  FindnumDe(0,0,nv,nvp,Decal,ErrDe)
              kw = 0
c----------------------------------------------------------
c  Calculate Morse parameters :
c----------------------------------------------------------
  20       al0 = Re * dsqrt( 0.50*amu*Wee*Wee/De )
         wexem = al0*al0/(2.0*amu*Re*Re)
c----------------------------------------------------------
c  Calculate vib-rot energies using vib-rot constants
c----------------------------------------------------------
           if (wei .eq. 0.0) wei = We
c----------------------------------------------------------
      if (kw .gt. 0) then
          write(38,21) w0,we0,wei,wex,cxe,wey,cye,wez,cze,
     #wet,cte,wes,cse,wer,cre
  21   format(/13x,'---  In subroutine  "calEvj"   --- ',/
     #/8x,'  w0 =',1PE18.10,
     #/8x,' we0 =',1PE18.10,
     #/8x,'  We =',1PE18.10,
     #/8x,'WeXe =',1PE18.10,5x,'cxe =',1pe10.2,
     #/8x,'WeYe =',1PE18.10,5x,'cye =',1pe10.2,
     #/8x,'WeZe =',1PE18.10,5x,'cze =',1pe10.2,
     #/8x,'WeTe =',1PE18.10,5x,'cte =',1pe10.2,
     #/8x,'WeSe =',1PE18.10,5x,'cse =',1pe10.2,
     #/8x,'WeRe =',1PE18.10,5x,'cre =',1pe10.2,/)
      endif 
c----------------------------------------------------------
c
         nv1 = nv+1
         nv2 = nv+50
c---
         if (nv .eq. 1 .and. nj(1) .eq. 1) then
           nv1 = 1
           nv2 = 1
         endif
c===
      do 25 i=1,nv2
        Ev(i) = 0.0d0
           bv = 1.0d0*( i - 1 )
c          bv = dfloat( i - 1 )
             if (nv2 .eq. 1) bv = v0
          bv0 = bv + 0.5d0
         bv02 =  bv0*bv0
         bv03 = bv02*bv0
         bv04 = bv03*bv0
         bv05 = bv04*bv0
         bv06 = bv05*bv0
         bv07 = bv06*bv0
c
         bv08 = bv07*bv0
         bv09 = bv08*bv0
         bv10 = bv09*bv0
         bv11 = bv10*bv0
         bv12 = bv11*bv0
         bv13 = bv12*bv0
         bv14 = bv13*bv0
         bv15 = bv14*bv0
c
C-
        if (nw0 .eq. 2) then
          if (w0 .lt. 0.0) then
            Ev2(i) =  cxe*w0 + (wei + we0)*bv0 
          else
            Ev2(i) =  w0 + (wei + we0)*bv0 
	    endif
        elseif (nw0 .eq. 1) then
          if (w0 .lt. 0.0) then
            Ev2(i) =  cxe*w0 + wei*bv0 
          else
            Ev2(i) =  w0 + wei*bv0 
	    endif
        elseif (nw0 .eq. 0) then
          Ev2(i) =  wei*bv0 
        endif
C-
         if (wex .lt. 0.0) then
           Ev2(i) =  Ev2(i) - cxe*wex*bv02 
         else
           Ev2(i) =  Ev2(i) - wex*bv02 
         endif
c
       IF (meny .eq. 0) THEN
           Ev3(i) =  Ev2(i) + aye*wey*bv03
           Ev4(i) =  Ev3(i) + aze*wez*bv04 
           Ev5(i) =  Ev4(i) + ate*wet*bv05 
           Ev6(i) =  Ev5(i) + ase*wes*bv06 
c
c--- Vibrational energies :
c
            Ev(i) =  Ev6(i) + are*wer*bv07
c-
c           Ev(i) =  Ev(i) + w08*bv08 + w09*bv09 + w10*bv10
c           Ev(i) =  Ev(i) + w11*bv11 + w12*bv12 + w13*bv13
c           Ev(i) =  Ev(i) + w14*bv14 + w15*bv15 
c---
       ELSE
           Ev3(i) =  Ev2(i) + aye*wey*bv03
           Ev4(i) =  Ev3(i) - aze*wez*bv04 
           Ev5(i) =  Ev4(i) + ate*wet*bv05 
           Ev6(i) =  Ev5(i) - ase*wes*bv06 
            Ev(i) =  Ev6(i) + are*wer*bv07
c-
c           Ev(i) =  Ev(i) - w08*bv08 + w09*bv09 - w10*bv10
c           Ev(i) =  Ev(i) + w11*bv11 - w12*bv12 + w13*bv13
c           Ev(i) =  Ev(i) - w14*bv14 + w15*bv15 
       ENDIF
c--------------------------------------
c Wee -- Original input We;
c  We -- Original input We if Mwe = 0 ;
c  We -- We + bryd  if Mwe > 0 .
c--- Approximate VIBrational energies :
c--------------------------------------
        if (wex .lt. 0.0) then
          Eu(i) = We*bv0 +  wex*bv02
        else
          Eu(i) = We*bv0 -  wex*bv02
        endif
        if (WeXe .lt. 0.0) then
          Ex(i) = Wee*bv0 + WeXe*bv02
        else
          Ex(i) = Wee*bv0 - WeXe*bv02
        endif
          Em(i) = Wee*bv0 - wexem*bv02
          Eh(i) = Wee*bv0 
c--- New terms of VIBrational energies :
        if (w0 .lt. 0.0) then
          Ew(i) = cxe*w0 + we0*bv0
        else
          Ew(i) = w0 + we0*bv0
        endif
c--- Ev's without NEW terms :
        Ew0(i) = Ev(i) - Ew(i)
        Eu0(i) = Ew0(i)
c----------------------------------------------------------
        if (nm .gt. 0 .or. kw .gt. 0) goto 25
c----------------------------------------------------------
            nj0 = nj(i)
        if (nj(1) .gt. 0 .and. i .le. nv1) then
          do j=1,nj0
              Evj(i,j) = 0.0d0
             bj = 1.0d0*( j - 1.0 )
               if (nv2 .eq. 1) bj = r0
            bj1 = bj*(bj + 1.0)
            bj0 = bj*(bj + 1.0) - lda*lda
            bj2 = bj0*bj0
            ejb = abe*Bee*bj0 -  aae*ale*bj0*bv0 
            ejb = ejb + age* gae*bj0*bv02 
            ejb = ejb - ae3*eta3*bj0*bv03 - ae4*eta4*bj0*bv04
            ejb = ejb - ae5*eta5*bj0*bv05 - ae6*eta6*bj0*bv06
            ej1 = ejb - ae7*eta7*bj0*bv07
            ej1a= ej1 - abe*Bee*bj0
c--
c           ejb = ejb + ae3*eta3*bj0*bv03 + ae4*eta4*bj0*bv04
c           ejb = ejb + ae5*eta5*bj0*bv05 + ae6*eta6*bj0*bv06
c           ej1 = ejb + ae7*eta7*bj0*bv07
c
c--- Herzberg (1953) :
C           ejc = - ade*Dee*bj2 - abt*bete*bj2*bv0
c--- A. L. G. Rees [ Proc. Phys. Soc. (London)59,998(1947) ] :
            ejc =   ade*Dee*bj2 - abt*bete*bj2*bv0
c
            ejc = ejc + ax2*xsi2*bj2*bv02 
            ejc = ejc + ax3*xsi3*bj2*bv03  + ax4*xsi4*bj2*bv04
            ejc = ejc + ax5*xsi5*bj2*bv05  + ax6*xsi6*bj2*bv06
            ej2 = ejc + ax7*xsi7*bj2*bv07 
            ej2a= ej2 - ade*Dee*bj2
c
c--- VIBrational-ROTational energies :
            Evj(i,j) = Ev(i) + ej1 + ej2
c--- ROTational energies :
            if (i .eq. 1) Ej(j) = abe*Bee*bj0 + ade*Dee*bj2
c--- ROTational energies using INPUT rot. constants :
            if (i .eq. 1) Eri(j) = abe*Be*bj0 + ade*Der*bj2
c--- VIBrational-ROTational COUpling energies :
            Ewj(i,j) = ej1a + ej2a
c--- Other energy forms :
            Euj(i,j) = Ev(i) + ej1 
            aj1(i,j) = ej1 
            aj2(i,j) = ej2 
c--- Rigid rotor energies at R=Re :
            if (i .eq. 1) Er(j) = Brigid*bj*(bj+1.0)
          enddo
        endif
c
  25  continue
          if (nv2 .eq. 1) goto 900 
c----------------------------------------------------------
          if (nm .gt. 0 .or. kw .gt. 0) goto 100
c----------------------------------------------------------
c  Write out vib-rot energies
c----------------------------------------------------------
        write( 6,500) 
        write(35,500) 
        write(36,500) 
            Emax = 0.0
          sumdif = 0.0
              kk = 0
      do 30 i=1,nvp
          kv = i - 1
          difv = Ev(i) - Ev(i-1)
            if (i .eq. 1) difv = Ev(i)
          if ( Ev(i) .lt. Ev(i-1) .and. kk .eq. 0) then
              kk = 1
            write( 6,*)
            write(35,*) 
            write(36,*) 
          endif
        write(35,510) kv, Ev(i), Ev(i)*aucm, difv, difv*aucm
        write( 6,510) kv, Ev(i), Ev(i)*aucm, difv, difv*aucm
        write(36,510) kv, Ev(i)
          if ( Ev(i) .gt. Ev(i-1) )   Emax = Ev(i)
          if (difv .gt. 0.0) sumdif = sumdif + difv
  30  continue
        write( 6,530) sumdif, Emax
        write(35,530) sumdif, Emax
c       write(35,530) sumdif, Emax-Ev(1)
c
        write( 6,535) 
        write(35,535) 
          kk = 0
      do i=1,nvp
            kv = i - 1
          difv = Ev(i) - Ev(i-1)
            if (i .eq. 1) difv = Ev(1)
          if ( Ev(i) .lt. Ev(i-1) .and. kk .eq. 0) then
              kk = 1
            write( 6,*)
            write(35,*) 
          endif
        write( 6,510) kv, Ev(i), Ew(i), difv
        write(35,510) kv, Ev(i), Ew(i), difv
      enddo
c
        write( 6,560) nw0 
        write(35,560) nw0
          kk = 0
      do i=1,nvp
            kv = i - 1
          difv = Ew0(i) - Ew0(i-1)
            if (i .eq. 1) difv = Ew0(1)
          if ( Ev(i) .lt. Ev(i-1) .and. kk .eq. 0) then
              kk = 1
            write( 6,*)
            write(35,*)
          endif
        write( 6,510) kv, Ev(i), Ew0(i), Ew(i), difv
        write(35,510) kv, Ev(i), Ew0(i), Ew(i), difv
      enddo
c
        write( 6,565)
        write(35,565)
          kk = 0
      do i=1,nvp
            kv = i - 1
          difv = ( Ew0(i) - Ew0(i-1) )*aucm
            if (i .eq. 1) difv = Ew0(1)*aucm
          if ( Ev(i) .lt. Ev(i-1) .and. kk .eq. 0) then
              kk = 1
            write( 6,*)
            write(35,*)
          endif
            evcm = Ev(i)*aucm
            ewcm = Ew(i)*aucm
            e0cm = Ew0(i)*aucm
        write( 6,510) kv, evcm, e0cm, ewcm, difv
        write(35,510) kv, evcm, e0cm, ewcm, difv
      enddo
c
        write( 6,570)
	  write(35,570)
          Emaxa = 0.0
          Emaxb = 0.0
          sumdifa = 0.0
          sumdifb = 0.0
          kk = 0
      do i=1,nvp
            kv = i - 1
          if ( Ev(i) .lt. Ev(i-1) .and. kk .eq. 0) then
              kk = 1
            write( 6,*)
            write(35,*)
          endif
              difa = Eu(i) - Eu(i-1)
              difb = Ex(i) - Ex(i-1)
            if (i .eq. 1) difa = Eu(i)
            if (i .eq. 1) difb = Ex(i)
          if (difa .gt. 0.0) sumdifa = sumdifa + difa
          if (difb .gt. 0.0) sumdifb = sumdifb + difb
            if ( Eu(i) .gt. Eu(i-1) ) Emaxa = Eu(i)
            if ( Ex(i) .gt. Ex(i-1) ) Emaxb = Ex(i)
        write(35,510) kv, Ev(i), Ew(i), Eu(i), Ex(i)
        write( 6,510) kv, Ev(i), Ew(i), Eu(i), Ex(i)
      enddo
        write( 6,540) sumdifa, sumdifb, Emaxa, Emaxb
        write(35,540) sumdifa, sumdifb, Emaxa, Emaxb
c
        write( 6,544)
        write(35,546)
          kk = 0
         difa = 0.0
         difb = 0.0
      do i=1,nvp
            kv = i - 1
          if ( Ev(i) .lt. Ev(i-1) .and. kk .eq. 0) then
              kk = 1
            write( 6,*)
            write(35,*)
          endif
            evcm = Ev(i)*aucm
            e0cm = Ew0(i)*aucm
            eucm = Eu(i)*aucm
            excm = Ex(i)*aucm
            difau = Eu(i) - Ex(i)
            difcm = difau*aucm
            difa = difa + abs(difau)
            difb = difb + abs(difcm)
        write( 6,510) kv, evcm, e0cm, eucm, excm
        write(35,510) kv, eucm, excm, difcm, difau
      enddo
        write( 6,550) difb/nv ,  difa/nv
        write(35,550) difb/nv ,  difa/nv
c
c================================================================
        IF (kin .eq. 0) GOTO 100 
C  40   IF (kin .eq. 0) GOTO 100 
c----------------------------------------------------------------
          call getEvDe(nv,kw,lda,Lj,Ev)
c--------------------------------------------------
C   Set VIB. spectrum constants which CONVERGES.
c--------------------------------------------------
        if (kw .gt. 0) then 
           if (nw0 .ge. 1) then
             wei = wvib0(2)
             wex = wvib0(3)
             wey = wvib0(4)
             wez = wvib0(5)
             wet = wvib0(6)
             wes = wvib0(7)
             wer = wvib0(8)
           else
             wei = wvib0(1)
             wex = wvib0(2)
             wey = wvib0(3)
             wez = wvib0(4)
             wet = wvib0(5)
             wes = wvib0(6)
             wer = wvib0(7)
C-
c            w08 = wvib0(8)
c            w09 = wvib0(9)
c            w10 = wvib0(10)
c            w11 = wvib0(11)
c            w12 = wvib0(12)
c            w13 = wvib0(13)
c            w14 = wvib0(14)
c            w15 = wvib0(15)
C-
           endif
             w0  = wvib0(30)
             we0 = wvib0(31)
          goto 20
	  endif
c================================================================
        write( 6,590)
        write(35,590)
          kk = 0
         esum= 0.0
C     do i=1,nv
      do i=1,nvd
            kv = i - 1
          if ( Ev(i) .lt. Ev(i-1) .and. kk .eq. 0) then
              kk = 1
            write( 6,*)
            write(35,*)
          endif
c
C           dvcm = ( Ev(i) - Ew0(i) )*aucm
C           ducm = ( Eu(i) - Ew0(i) )*aucm
C           dxcm = ( Ex(i) - Ew0(i) )*aucm
C           dvcm = ( ( Ev(i)-Ev(i-1) ) - Ew0(i) )*aucm
c 
            e1(i)= ( Eu(i)-Eu(i-1) )*aucm
            e2(i)= ( ( Eu(i)-Eu(i-1) ) - Ew0(i) )*aucm
            e3(i)= 100.0 * abs( e2(i) )/( Ew0(i)*aucm )
            dxcm = ( ( Ex(i)-Ex(i-1) ) - Ew0(i) )*aucm
              if (i .eq. 1) dvcm = ( Ev(1) - Ew0(1) )*aucm
              if (i .eq. 1) e2(i)= ( Eu(1) - Ew0(1) )*aucm
              if (i .eq. 1) dxcm = ( Ex(1) - Ew0(1) )*aucm
            if ( i .le. nvd ) esum = esum + e3(i)
        write( 6,510) kv, Ew0(i)*aucm, e1(i), e2(i), dxcm
        write(35,510) kv, Ew0(i)*aucm, e1(i), e2(i), dxcm
C       write(35,510) kv, Ew0(i)*aucm, dvcm,  e2(i), dxcm
      enddo
c 
        write( 6,594)
        write(35,594)
      do i=1,nvd
            kv = i - 1
        write( 6,510) kv, Ew0(i)*aucm, e1(i), e2(i), e3(i)
        write(35,510) kv, Ew0(i)*aucm, e1(i), e2(i), e3(i)
      enddo
        write( 6,596) esum/nvd
        write(35,596) esum/nvd
c
        write( 6,600)
        write(35,600)
          esum = 0.0
      do i=1,nvd
              kv = i - 1
            e1(i)= ( Ev(i)-Ev(i-1) )*aucm
            e2(i)= ( ( Ev(i)-Ev(i-1) ) - Ew0(i) )*aucm
            e3(i)= 100.0 * abs( e2(i) )/( Ew0(i)*aucm )
              if ( i .le. nvd ) esum = esum + e3(i)
        write( 6,510) kv, Ew0(i)*aucm, e1(i), e2(i), e3(i)
        write(35,510) kv, Ew0(i)*aucm, e1(i), e2(i), e3(i)
      enddo
        write( 6,602) esum/nvd
        write(35,602) esum/nvd
c
        write( 6,604)
        write(35,604)
          esum = 0.0
      do i=1,nvd
              kv = i - 1
            e3(i)= 100.0 * abs( Ew(i) - Ev(i) )/Ew(i)
              if ( i .le. nvd ) esum = esum + e3(i)
        write( 6,510) kv, Ew(i)*aucm, Ev(i)*aucm, e3(i)
        write(35,510) kv, Ew(i)*aucm, Ev(i)*aucm, e3(i)
      enddo
        write( 6,606) esum/nvd
        write(35,606) esum/nvd
c
        write( 6,634)
        write(35,634)
          esum = 0.0
      do i=1,nvd
              kv = i - 1
            e3(i)= 100.0 * abs( Ew(i) - Eu0(i) )/Ew(i)
              if ( i .le. nvd ) esum = esum + e3(i)
        write( 6,510) kv, Ew(i)*aucm, Eu0(i)*aucm, e3(i)
        write(35,510) kv, Ew(i)*aucm, Eu0(i)*aucm, e3(i)
      enddo
        write( 6,636) esum/nvd
        write(35,636) esum/nvd
c===
 100  if (wex .lt. 0.0) then
        write( 6,620) cre*wex,WeXe,wexem,
     #wex*aucm,WeXe*aucm,wexem*aucm 
        write(35,620) cre*wex,WeXe,wexem,
     #wex*aucm,WeXe*aucm,wexem*aucm 
      else
        write( 6,620) wex,WeXe,wexem,
     #wex*aucm,WeXe*aucm,wexem*aucm 
        write(35,620) wex,WeXe,wexem,
     #wex*aucm,WeXe*aucm,wexem*aucm 
	endif
         kk = 0
      do i=1,nv2
         kv = i - 1
        ev0 = Ev(i)
        ex0 = Ex(i)
        em0 = Em(i)
        eh0 = Eh(i)
          if ( Ev(i) .lt. Ev(i-1) ) ev0 = 0.0
          if ( Ex(i) .lt. Ex(i-1) ) ex0 = 0.0
          if ( Em(i) .lt. Em(i-1) ) em0 = 0.0
c         if ( Eh(i) .lt. Eh(i-1) ) eh0 = 0.0
          if ( Ev(i) .lt. Ev(i-1) .and. kk .eq. 0) then
              kk = 1
            write( 6,*)
            write(35,*) 
          endif
        write(35,510) kv, ev0, ex0, em0, eh0
        write( 6,510) kv, ev0, ex0, em0, eh0
      enddo
C
        write( 6,610)
        write(35,610)
        write(38,642)
      do k=1,6
            if (k .eq. 1) write( 6,611) 
            if (k .eq. 1) write(35,611) 
            if (k .eq. 2) write( 6,612) 
            if (k .eq. 2) write(35,612) 
            if (k .eq. 2) write(38,644) 
            if (k .eq. 3) write( 6,613) 
            if (k .eq. 3) write(35,613) 
            if (k .eq. 3) write(38,645) 
            if (k .eq. 4) write( 6,614) 
            if (k .eq. 4) write(35,614) 
            if (k .eq. 4) write(38,646) 
            if (k .eq. 5) write( 6,615) 
            if (k .eq. 5) write(35,615) 
            if (k .eq. 5) write(38,647) 
            if (k .eq. 6) write(38,648) 
              esum = 0.0
	  do i=1,nvp
            kv = i-1
          if (k .eq. 1) then
              error = 100.0*( Ev3(i) - Ev2(i) )/Ev2(i)
            write( 6,510) kv, Ev2(i)*aucm, Ev3(i)*aucm, error
            write(35,510) kv, Ev2(i)*aucm, Ev3(i)*aucm, error
C           write(38,510) kv, Ev2(i)*aucm, Ev3(i)*aucm, error
            write(38,510) kv, Ev2(i)
	    elseif(k .eq. 2) then
              error = 100.0*( Ev4(i) - Ev3(i) )/Ev3(i)
            write( 6,510) kv, Ev3(i)*aucm, Ev4(i)*aucm, error
            write(35,510) kv, Ev3(i)*aucm, Ev4(i)*aucm, error
C           write(38,510) kv, Ev3(i)*aucm, Ev4(i)*aucm, error
            write(38,510) kv, Ev3(i)
	    elseif(k .eq. 3) then
              error = 100.0*( Ev5(i) - Ev4(i) )/Ev4(i)
            write( 6,510) kv, Ev4(i)*aucm, Ev5(i)*aucm, error
            write(35,510) kv, Ev4(i)*aucm, Ev5(i)*aucm, error
C           write(38,510) kv, Ev4(i)*aucm, Ev5(i)*aucm, error
            write(38,510) kv, Ev4(i)
	    elseif(k .eq. 4) then
              error = 100.0*( Ev6(i) - Ev5(i) )/Ev5(i)
            write( 6,510) kv, Ev5(i)*aucm, Ev6(i)*aucm, error
            write(35,510) kv, Ev5(i)*aucm, Ev6(i)*aucm, error
C           write(38,510) kv, Ev5(i)*aucm, Ev6(i)*aucm, error
            write(38,510) kv, Ev5(i)
	    elseif(k .eq. 5) then
              error = 100.0*abs(  Ev(i) - Ev6(i) )/Ev6(i)
            write( 6,510) kv, Ev6(i)*aucm, Ev(i)*aucm, error
            write(35,510) kv, Ev6(i)*aucm, Ev(i)*aucm, error
C           write(38,510) kv, Ev6(i)*aucm, Ev(i)*aucm, error
            write(38,510) kv, Ev6(i)
	    elseif(k .eq. 6) then
            write(38,510) kv, Ev(i)
	    endif
            if (k .lt. 6) esum = esum + error
        enddo
            if (k .lt. 6) then
		  write( 6,618) esum/nvp
		  write(35,618) esum/nvp
            endif
      enddo
C-
        IF (mtyp .eq. 3) GOTO 900
        IF ( nm .gt. 0 ) GOTO 900
c-
c====================================================
c-- bjj is classicle rotational ANGular momentum
c-- bvj is rotational ANGular velocity
c
        write( 6,710) rinert, Brigid
        write(35,710) rinert, Brigid
        write(36,710) rinert, Brigid
      do j=1,nj(1)
           bj = j*1.0 - 1.0
          bjj = dsqrt( bj*(bj + 1.0) )
          bvj = bjj/rinert
        write( 6,510) j-1, Ej(j), Er(j), bjj, bvj
        write(35,510) j-1, Ej(j), Er(j), bjj, bvj
        write(36,510) j-1, Ej(j)
      enddo
c
        write( 6,720)
        write(35,720)
	    sumj = 0.0
      do j=1,nj(1)
        if (j .eq. 1) then
	    bjdif = 0.0
        else
          bjdif = 100.0 * abs( Ej(j) -  Er(j) )/Ej(j)
        endif
	        sumj = sumj + bjdif
        write( 6,510) j-1, Ej(j), Er(j), bjdif
        write(35,510) j-1, Ej(j), Er(j), bjdif
      enddo
        write( 6,730) sumj/nj(1)
        write(35,730) sumj/nj(1)
c
        write( 6,740)
        write(35,740)
	    sumj = 0.0
      do j=1,nj(1)
        if (j .eq. 1) then
	    bjdif = 0.0
        else
          bjdif = 100.0 * abs( Ej(j) -  Er(j) )/Ej(j)
        endif
	        sumj = sumj + bjdif
        write( 6,510) j-1, Ej(j), Eri(j), bjdif
        write(35,510) j-1, Ej(j), Eri(j), bjdif
      enddo
        write( 6,730) sumj/nj(1)
        write(35,730) sumj/nj(1)
c
      if (nj(1) .gt. 0) then
          write( 6,625) ade 
          write(35,625) ade
          write(36,625) ade
            sumdif = 0.0
             difj0 = Evj(1,1)
            Evjmax = 0.0
            jj = 0
        do i=1,nv
            kv = i - 1
            nj0 = nj(i)
             kk = 0
		   if (Evj(i,1) .le. Emax) then
		     if (Evj(i,1) .gt. Evjmax) then
		       Evjmax = Evj(i,1)
                 endif
               endif
          do j=1,nj0
            kj = j - 1
              if (j .ge. 2) difvj = Evj(i,j) - Evj(i,j-1)
              if (j .eq. 1) difvj = difj0
             if ( Evj(i,j) .lt. Evj(i+1,1) ) then
		   sumdif = sumdif + difvj
             endif
             if ( Evj(i,j) .gt. Evj(i+1,1) .and. kk .eq. 0 ) then
                   difj0 = Evj(i+1,1) - Evj(i,j-1)
 	         if (Evj(i,j-1) .gt. Evjmax) then
		     if (Evj(i,j-1) .le. Emax)  Evjmax = Evj(i,j-1)
		   endif
                   kk = 1
                 write( 6,*) 
                 write(35,*) 
                 write(36,*) 
             endif
           write(35,650) kv,kj,Evj(i,j),Evj(i,j)*aucm,difvj,difvj*aucm
           write( 6,650) kv,kj,Evj(i,j),Evj(i,j)*aucm,difvj,difvj*aucm
           write(36,650) kv,kj,Evj(i,j)
             if ( Evj(i,j) .eq. Emax)  jj = 1
             if ( Evj(i,j) .eq. Emax)  sumdif = sumdif + difvj
          enddo
	      if (jj .eq. 1) then
              write( 6,*) 
              write(35,*) 
              write(36,*) 
		else
              write( 6,630)
              write(35,630)
              write(36,630)
		endif
              if (i .eq. nv)  sumdif = sumdif - (sumdif - Evjmax)
        enddo
	      dif = sumdif - Evjmax
          write(35,640) sumdif, Evjmax, dif
          write( 6,640) sumdif, Evjmax, dif
c
          write(35,660) 
          write( 6,665) 
            jj = 0
        do i=1,nv
            kv = i - 1
            nj0 = nj(i)
             kk = 0
          do j=1,nj0
              kj = j - 1
            if ( Evj(i,j) .gt. Evj(i+1,1) .and. kk .eq. 0) then
                kk = 1
              write(35,*)
              write( 6,*)
            endif
            write(35,650) kv,kj,Evj(i,j),Ev(i),aj1(i,j)+aj2(i,j)
            write( 6,650) kv,kj,Evj(i,j),Ewj(i,j),Ev(i),Ej(j)
              if ( Evj(i,j) .eq. Emax) jj = 1
          enddo
            if (jj .eq. 1) then
              write( 6,*)
              write(35,*)
            else
              write( 6,630)
              write(35,630)
            endif
        enddo
c
          write(35,670) 
        do i=1,nv
            kv = i - 1
            nj0 = nj(i)
          do j=1,nj0
            kj = j - 1
            write(35,650) kv,kj,Evj(i,j),Ev(i),aj1(i,j),aj2(i,j)
          enddo
            write(35,*) 
        enddo
c
          write(35,680) 
        do i=1,nv
            kv = i - 1
            nj0 = nj(i)
          do j=1,nj0
            kj = j - 1
            write(35,650) kv,kj,Evj(i,j),Euj(i,j),aj1(i,j)
          enddo
            write(35,*) 
        enddo
c
      endif
c
c----------------------------------------------------------
 400  format(///6x,'--- Ro-vibrational constants & energies --- ',/
     #/2x,'All calculated constants are the functions of the input ', 
     #/2x,'"We", "Re", & "mu",  and are based on the perturbation ',
     #/2x,'theory,  except that "We" is the input value.',
     #///8x,'**  The vibrational constants (in a.u.) are  ** ',
     #/6x,'evaluated using force constants from ANAlytical deriv.', 
C    #/6x,'evaluated using force constants from NUMerical deriv.', 
     #//6x,'  { If mtyp = 1, Consts(calc.) == Consts(input) } ',
     #//10x,'[ Error% = 100.0*(|Calc.| - |Input|)/Calc. ] ')
 410  format(/5x,'          Calculated              Input ',
     #13x,'Error% ',
     #//3x,' W0  = ', 1PE20.12,/3x,' We0 = ', 1PE20.12,
     # /3x,' We  = ', 1PE20.12,2x,1PE20.12,
     # /3x,'WeXe = ', 1PE20.12,2x,1PE20.12,2x,1PE13.5,
     # /3x,'WeYe = ', 1PE20.12,2x,1PE20.12,2x,1PE13.5,
     # /3x,'WeZe = ', 1PE20.12,2x,1PE20.12,2x,1PE13.5,
     # /3x,'WeTe = ', 1PE20.12,2x,1PE20.12,2x,1PE13.5,
     # /3x,'WeSe = ', 1PE20.12,2x,1PE20.12,2x,1PE13.5,
     # /3x,'WeRe = ', 1PE20.12,2x,1PE20.12,2x,1PE13.5 )
 420  format(//13x,'The vibrational constants (in cm-1) are : ')
 430  format(///12x,'## The rotational constants (in a.u.) are : ## ')
 435  format(///12x,'Rotational consts (in a.u.) from f_(n; cal.) : ')
 440  format(/5x,'             Calculated               Input ',
     #13x,'Error% ',
     #//3x,'     Be = ', 1PE20.12,2x,1PE20.12,2x,1PE13.5,
     # /3x,'Alpha_e = ', 1PE20.12,2x,1PE20.12,2x,1PE13.5,
     # /3x,'Gamma_e = ', 1PE20.12,2x,1PE20.12,2x,1PE13.5,
     # /3x,' Eta_e3 = ', 1PE20.12, /3x,' Eta_e4 = ', 1PE20.12,
     # /3x,' Eta_e5 = ', 1PE20.12, /3x,' Eta_e6 = ', 1PE20.12,
     # /3x,' Eta_e7 = ', 1PE20.12,
     #//3x,'    Dee = ', 1PE20.12,2x,1PE20.12,2x,1PE13.5, 
     # /3x,' Beta_e = ', 1PE20.12,2x,1PE20.12,2x,1PE13.5,
     # /3x,' Xsi_e2 = ', 1PE20.12, /3x,' Xsi_e3 = ', 1PE20.12,
     # /3x,' Xsi_e4 = ', 1PE20.12, /3x,' Xsi_e5 = ', 1PE20.12,
     # /3x,' Xsi_e6 = ', 1PE20.12, /3x,' Xsi_e7 = ', 1PE20.12 )
 450  format(//15x,'The rotational constants (in cm-1) are : ')
 460  format(///15x,'    Check VIBRATIONAL constants : ',
     #//8x,"[ As V(R)=De, De can be calculated by Rees's Eq. ",
     # /8x,"            from  VIBrational constants.         ",/
     #//8x,"        *** Quadratic potential form : ***  ",/
     # /8x,"  De_2c -- De from calculated vibrational constants; ",
     # /8x,"  De_2i -- De from  inputted  vibrational constants. ",
     #//8x,"                De_2 = We*We/(4*Wexe)              ",/
     #//8x,"        ***  CUBIC   potential form :  ***  ",/
     # /8x,"  De_3c -- De from calculated vibrational constants; ",
     # /8x,"  De_3i -- De from  inputted  vibrational constants. ",
     #//8x,"     De3(a) = 2*[dsqrt(Wexe**2 - 3*We*Weye)]**3      ",
     # /8x,"     De3(b) =  Wexe*(2*Wexe**2 - 9*We*Weye)          ",
     #//8x,"       De_3 = [De3(a) - De3(b)]/(27*Weye**2)         ",
     #//8x,"     IF  term = Wexe**2 - 3*We*Weye  < 0.0           ",
     # /8x,"       dsqrt( term ) == complex number !             ",
     #//8x,"       De_2   =       De_2c  OR  De_2i               ",
     # /8x,"       De_3   =       De_3c  OR  De_3i               ",
     # /8x,"       De3(a) =       De_3ca OR  De_3ia              ",
     # /8x,"       De3(b) =       De_3cb OR  De_3ib              ",
     #//8x,"    Set De_3i & D3i% = 0.0 if inputted WeYe is 0.0 . ",
     #///8x,"        *** Quadrupole potential form : *** ",/
     # /8x,"  De_4c -- De from calculated vibrational constants; ",
     # /8x,"  De_4i -- De from  inputted  vibrational constants. ",
     #//8x,"      delc = 9.0*weye**2 + 24.0*weze*wexe  >= 0  ",
     # /8x,"      deli = 9.0*WeYe**2 + 24.0*WeZe*WeXe  >= 0  ",
     #//8x,"    yc = ( -3.0*weye +/- dsqrt(delc) )/(12.0*weze) ",
     # /8x,"    yi = ( -3.0*WeYe +/- dsqrt(deli) )/(12.0*WeZe) ",
     #///8x,"        ***  Quintic  potential form : *** ",/
     # /8x,"  De_5c -- De from calculated vibrational constants; ",
     # /8x,"  De_5i -- De from  inputted  vibrational constants. ",
     #//8x,"  Higher_power potential form is much more accurate. ",
     #//8x,"       Dei = Input 'TRUE' dissociation energy De :   ",
     #//8x,"           D2c% = 100 * | De_2c - Dei |/Dei ;   ",
     # /8x,"           D2i% = 100 * | De_2i - Dei |/Dei ;   ",
     # /8x,"        D3c% = 100 * | real(De_3c - Dei) |/Dei ;   ",
     # /8x,"        D3i% = 100 * | real(De_3i - Dei) |/Dei ;   ",
     # /8x,"        D4c% = 100 * | real(De_4c - Dei) |/Dei ;   ",
     # /8x,"        D4i% = 100 * | real(De_4i - Dei) |/Dei .   ",
     # /8x,"        D5c% = 100 * | real(De_5c - Dei) |/Dei ;   ",
     # /8x,"        D5i% = 100 * | real(De_5i - Dei) |/Dei .       ]",
     #//18x,"    Dei == De = ",f12.8,'  a.u.',
     #//18x,"    From ANALYTICAL calculations : ",   
     #//8x,"n    De_nc(a.u)",6x,"D_nc%",9x,"De_ni(a.u)",6x,"D_ni%",/
     #4(/7x,i2,4(3x,1pe12.6),) )
 470  format(///15x,'    Check CUBIC form calculations : ',
     #//12x," Using the calculated constants Wexe and Weye ",
     #//11x,"       Wexe**2(a.u.)     -3*We*Weye(a.u.)   ",
     #/17x,14(1H-),5x,16(1H-)
     #/15x,2(1PE16.8,4x),
     #///4x,"De_3cb(a.u)",9x,"De_3ca(a.u)",16x,"De_3c(a.u)",
     #/3x,13(1H-),2x,24(1H-),2x,25(1H-)
     #/7x,"Real",10x,"Real",7x,"Imaginary",6x,"Real",8x,"Imaginary",
     #/2x,1PE14.6,2(1PE13.5),1x,2(1PE13.5) )
 475  format(///15x,'    Check CUBIC form calculations : ',
     #//15x," Using the INPUT vibrational constants ",
     #//11x,"       WeXe**2(a.u.)     -3*We*WeYe(a.u.)   ",
     #/17x,14(1H-),5x,16(1H-)
     #/15x,2(1PE16.8,4x),
     #///4x,"De_3ib(a.u)",9x,"De_3ia(a.u)",16x,"De_3i(a.u)",
     #/3x,13(1H-),2x,24(1H-),2x,25(1H-)
     #/7x,"Real",10x,"Real",7x,"Imaginary",6x,"Real",8x,"Imaginary",
     #/2x,1PE14.6,2(1PE13.5),1x,2(1PE13.5) )
 480  format(///12x,'    Check QUADRUPOLE form calculations : ',
     #/24x,"  For  delc  <  0.0 ",
     #//12x," Using calculated constants Wexe, Weye & Weze ",
     #//15x,"delc(a.u.)",3x,"Root_yc(a.u.)",2x,"De_4c(a.u)",
     #/13x,13(1H-),2x,12(1H-),2x,12(1H-)
     #/18x,"Real",10x,"Real",10x,"Real",
     #/12x,1PE14.6,1x,1PE13.5,1x,1PE13.5 )
 485  format(///12x,'    Check QUADRUPOLE form calculations : ',
     #/24x,"  For  deli  <  0.0 ",
     #//12x," Using INPUT constants WeXe, WeYe & WeZe ",
     #//15x,"deli(a.u.)",3x,"Root_yi(a.u.)",2x,"De_4i(a.u)",
     #/13x,13(1H-),2x,13(1H-),2x,12(1H-)
     #/18x,"Real",10x,"Real",10x,"Real",
     #/12x,1PE14.6,1x,1PE13.5,1x,1PE13.5 )
 500  format(///17x,'=== The VIBrational energies are : === ',
     #//19x,'    [  Ev_dif = E(v) - E(v-1)  ] ',
     #//5x,'v      E(v; a.u.)        E(v; cm-1)     Ev_dif(a.u)',
     #'   Ev_dif(cm-1)',/)
 510  format(3x,i3,1PE18.9,1PE18.9,x,1PE13.5,x,1PE13.5)
 520  format(1PE18.9)
 530  format(/12x,'The SUM of  Ev_dif(a.u) = ',1PE16.8,
     #/12x,'Maximum energy  Ev(max) = ',1PE16.8,/)
C    #/12x,'Ev(maximum)  -  Ev(v=0) = ',1PE16.8,/)
 535  format(//5x,'v      E(v; a.u.)      Enew(v; a.u.)',
     #'    Ev_dif(a.u.)',/)
 540  format(/2x,'SUM of Ea_dif(a.u) =',1PE13.6,
     # 3x,'SUM of Eb_dif(a.u) =',1PE13.6,
     #/2x,'Maxi. Ea(v)=Ea(max)=',1PE13.6,
     # 3x,'Maxi. Eb(v)=Eb(max)=',1PE13.6,/)
 544  format(//5x,"v     E(v; cm-1)       E'(v; cm-1)",
     #"      Ea(v; cm-1)    Eb(cm-1)",/)
 546  format(//5x,"v     Ea(v; cm-1)       Eb(v; cm-1)",
     #"     Ea-Eb(cm-1)   Ea-Eb(a.u.)",/)
 550  format(/15x,'Average of |Ea(cm-1) - Eb(cm-1)| =',1PE13.6,
     #/15x,'Average of |Ea(a.u.) - Eb(a.u.)| =',1PE13.6,/)
 560  format(///8x,"[  E'(v) = E(v) - Enew(v);  ",
     #"E(v) =/= E'(v)  for nw0 = 1 ; ",
     #/12x,"E(v) = E(v) - Enew(v);  E(v) === E'(v)  for nw0 = 0  ]",
     #/31x,"nw0 =",i2,/
     #/5x,"v     E(v; a.u.)        E'(v; a.u.)    Enew(v; a.u) ",
     #"  E'v_dif(au)",/)
 565  format(//5x,"v     E(v; cm-1)        E'(v; cm-1)",
     #"    Enew(v;cm-1)  E'v_dif(cm-1)",/)
 570  format(//16x,'Comparing vibrational energies (I) : ',
     #//12x,'[ Enew(v) = w0 + we0*(v + 0.5) ; ',
     # /12x,'    Ea(v) = We*(v + 0.5) - wexe*(v + 0.5)**2 ; ',
     # /12x,'    Eb(v) = We*(v + 0.5) - WeXe*(v + 0.5)**2 ; ',
     # /12x,'     E(v) = Enew(v) + Ea(v) - ...            ; ',
     #//12x,'          wexe is calculated by code; ',
     # /12x,'          WeXe is  the  INPUT  value;          ',
     # /12x,'       w0 & we0 are calculated NEW terms.        ',
     #//12x,'     If mtyp = 1, w0 & we0 are NOT defined.     ] ',
     #//5x,'v      E(v; a.u.)     Enew(v; a.u.)     Ea(v; a.u.)',
     #'   Eb(v; a.u.)',/)
 574  format(//8x,"Einp(v) == INPUT energies Ev's ; ",
     #/7x,"Edif_1(v) =   Einp(v)  -  Einp(v-1) ",
     #/7x,"Edif_2(v) =  Edif_1(v) - Edif_1(v-1) ",
     #//5x,'v     Einp(v; a.u.)   Edif_1(v; a.u.)  Edif_2(v; a.u.)',/)
 576  format(//5x,'v     Einp(v; cm-1)   Edif_1(v; cm-1) ',
     #' Edif_2(v; cm-1) ',/)
 580  format(///8x,'===== Following are generated from',
     #" readed Ev's =====",)
 582  format(/6x,'===== Finish calculations for',
     #" READED Ev's =====",/)
 584  format(///6x,'===== Following are generated from',
     #" CALCULATED Ev's =====",)
 586  format(/6x,'===== Finish calculations for',
     #" CALCULATED Ev's =====",)
C590  format(///25x,"  Einp(v) is the INPUT Ev's ",
C    #/27x,"Edifv(v) =  E(v) - Einp(v) ",
C    #/27x,"Edifa(v) = Ea(v) - Einp(v) ",
C    #/27x,"Edifb(v) = Eb(v) - Einp(v) ",/
C    #/5x,"v    Einp(v; cm-1)    Edifv(v; cm-1)",
C    #"    Edifa(cm-1)   Edifb(cm-1)",/)
 590  format(///18x,"DEin(v) is the INPUT differencies Edif_i(v) ",
C    #/18x,"Edelv(v) =  [  E(v) -  E(v-1) ] - DEin(v) ",
     #/18x,"Edifa(v) =  [ Ea(v) - Ea(v-1) ] ",
     #/18x,"Edela(v) =  [ Ea(v) - Ea(v-1) ] - DEin(v) ",
     #/18x,"Edelb(v) =  [ Eb(v) - Eb(v-1) ] - DEin(v) ",/
C    #/5x,"v    Einp(v; cm-1)    Edelv(v; cm-1)",
     #/5x,"v   Edif_i(v; cm-1)   Edifa(v; cm-1)",
     #"    Edela(cm-1)   Edelb(cm-1)",/)
 594  format(///17x,"ERRORa(v)% = 100.0*| Edela(v) |/Edif_i(v) : ",/
     #/5x,"v   Edif_i(v; cm-1)   Edifa(v; cm-1)",
     #"    Edela(cm-1)   ERRORa(cm-1)",/)
 596  format(/5x,'Average error for QUADRATIC Ev,',
     #'  ERRORa(v)% = ',1PE16.8,'%',/)
 600  format(//18x,"DEin(v) is the INPUT differencies Edif_i(v) ",/
     #/18x," Edifv(v) =  [ E(v) - E(v-1) ] ",
     #/18x," Edelv(v) =  [ E(v) - E(v-1) ] - DEin(v) ",
     #/18x,"ERROR(v)% = 100.0 * | Edelv(v) |/Edif_i(v) : ",/
     #/5x,"v   Edif_i(v; cm-1)   Edifv(v; cm-1)",
     #"    Edelv(cm-1)   ERROR(cm-1)",/)
 602  format(/5x,'Average error for FULL Ev difference,',
     #' ERROR(v)% =',1PE14.6,'%',/)
 604  format(//12x,"Einp(v) is the INPUT vib. energies ",/
     #/8x,"ERROR(v)% = 100.0*| Einp(v)-Ev(v)|/Einp(v) : ",/
     #/5x,"v    Einp(v; cm-1)      Ev(v; cm-1)     ERROR(cm-1) ",/)
 606  format(/4x,'Average error for FULL energies Ev,',
     #' ERROR(v)% =',1PE11.4,'%',/)
 610  format(///11x,'Comparing vibrational energies (III) : ',/
     #/12x,"Edifv(i) =  Ev_(n)[i] - Ev_(n-1)[i] ",
     #/8x,'Error_v(i)% = 100.0*{ Edifv(i) }/Ev_(n-1)[i] : ',)
 611  format(//11x,'Ev2(i) => Energy expanded to Wexe term,',
     #/11x,'Ev3(i) => Energy expanded to Weye term.',/
     #/5x,"v      Ev2(cm-1)         Ev3(cm-1)       Error_v%",/)
 612  format(//11x,'Ev3(i) => Energy expanded to Weye term,',
     #/11x,'Ev4(i) => Energy expanded to Weze term.',/
     #/5x,"v      Ev3(cm-1)         Ev4(cm-1)       Error_v%",/)
 613  format(//11x,'Ev4(i) => Energy expanded to Weze term,',
     #/11x,'Ev5(i) => Energy expanded to Wete term.',/
     #/5x,"v      Ev4(cm-1)         Ev5(cm-1)       Error_v%",/)
 614  format(//11x,'Ev5(i) => Energy expanded to Wete term,',
     #/11x,'Ev6(i) => Energy expanded to Wese term.',/
     #/5x,"v      Ev5(cm-1)         Ev6(cm-1)       Error_v%",/)
 615  format(//11x,'Ev6(i) => Energy expanded to Wese term,',
     #/11x,'Ev7(i) => Energy expanded to Were term.',/
     #/5x,"v      Ev6(cm-1)         Ev7(cm-1)       Error_v%",/)
 618  format(/7x,'Average percent error  ERROR(v)% =',1PE11.4,'%',/)
 620  format(//17x,'Comparing vibrational energies (II) : ',
     #//12x,' [  Eb(v) = We*(v + 0.5) - WeXe *(v + 0.5)**2 ; ',
     # /12x,'    Em(v) = We*(v + 0.5) - wexem*(v + 0.5)**2 ; ',
     # /12x,'    Eh(v) = We*(v + 0.5)                      . ',
     #//12x,'          WeXe   is  the   INPUT  value;        ',
     # /12x,'          wexem  is  from  Morse formulae.       ',
     #//12x,"    Eb is the approximate exp't vib. energies; ",
     # /12x,"    Em is the  Morse   vibrational   energies; ",
     # /12x,"    Eh is the   SHO    vibrational   energies.   ] ",
     #//10x,"        Calculated        Input         Morse   ",    
     # /10x,"           wex            WeXe          wexem   ",
     # /9x,"WeXe =",3(1PE14.6,x),'  a.u. ',
     # /15x,3(1PE14.6,x),'  cm-1 ',
     #//5x,'v      E(v; a.u.)      Eb(v; a.u.)      Em(v; a.u.)',
     #'   Eh(v; a.u.)',/)
 625  format(///21x,'The VIB-ROTational energies are : ',
     #//15x,'[ A. L. G. Rees, Proc. Phys. Soc.(London) ',
     # /15x,'             59, 998(1947) :              ',
     # /15x,'  E(v,j) = ..., + D_e*J*J*(J+1)**2 + ... ',
     #//15x,'       G. Herzberg, (1953) :              ',
     # /15x,'  E(v,j) = ..., - D_e*J*J*(J+1)**2 + ... ',
     #//16x,'       Evj_dif = E(v,j) - E(v,j-1)         ] ',
     #//15x,"{ ade =  1.0, use   Rees's   definition; ",
     # /15x,"      = -1.0, use Herzberg's definition.    } ",
     #//31x,"  ade =",f5.1,//
     #2x,'v  j    E(v,j; a.u.)      E(v,j; cm-1)    Evj_dif(a.u)',
     #'  Evj_dif(cm-1)',/)
 630  format(13x,'-----',13x,'-----',12x,'-----',9x,'-----',/)
 634  format(//17x,"Eu0(v) = Ev(v) - Enew(v) ",/
     #/8x,"ERROR(v)% = 100.0*| Einp(v)-Eu0(v)|/Einp(v) : ",/
     #/5x,"v    Einp(v; cm-1)      Eu0(v;cm-1)     ERROR(cm-1) ",/)
 636  format(/4x,'Average error for FULL energies Eu0,',
     #' ERROR(v)% =',1PE11.4,'%',/)
 640  format(/12x,'The SUM of Evj_dif(a.u) = ',1PE16.8,
     #/12x,'Maximum energy Evj(max) = ',1PE16.8,/
     #/12x,'Evj_dif(a.u) - Evj(max) = ',1PE16.8,//
     #/6x,'If  Evj_dif(a.u) =/= Evj(max),  Check  nj  in each ',
     #'v state !',
     #/11x,'Set  nj >= 1 + j_max (of vib. state) = 1 + (nj-1) '/)
 642  format(//1x,"Vibrational energies Ev's in various form : ",/
     #/1x,'(Generated from calculated vib. constants)',/
     #/2x,'Ev2(i) => Energy expanded to Wexe term;',
     #/2x,'Ev3(i) => Energy expanded to Weye term;',
     #/2x,'Ev4(i) => Energy expanded to Weze term;',
     #/2x,'Ev5(i) => Energy expanded to Wete term;',
     #/2x,'Ev6(i) => Energy expanded to Wese term;',
     #/2x,'Ev7(i) => Energy expanded to Were term.',/
     #/5x,"v      Ev2(a.u.)   ",/)
 644  format(//5x,"v      Ev3(a.u.)   ",/)
 645  format(//5x,"v      Ev4(a.u.)   ",/)
 646  format(//5x,"v      Ev5(a.u.)   ",/)
 647  format(//5x,"v      Ev6(a.u.)   ",/)
 648  format(//5x,"v      Ev7(a.u.)   ",/)
 650  format(2i3,1PE18.9,1PE18.9,x,1PE13.5,x,1PE13.5)
 660  format(///13x,'Contributions to E_vj are (Part I) : ',
     #//7x,'[   E(v,j) = E(v) + Ej ;     Ej = Ej1 + Ej2  ; ',
     # /7x,'       Ej1 = E{v;   j*(j+1)-lda*lda } ;    ',/
     # /7x,'       Ej2 = E{v; ( j*(j+1)-lda*lda )**2 }     ] ',/
     #/2x,'v  j    E(v,j; a.u.)        E(v;a.u.)       Ej(a.u.) ',/)
 665  format(///18x,'Separate contributions to E_vj  : ',
     #//12x,'[   E(v,j) = E(v-j_coupling) + E(v) + E(j) ; ',
     #//12x,'    E(v-j_coup) = Vib-Rot_coupling energies; ',
     # /12x,'           E(v) =    Vibrational   energies; ',
     # /12x,'           E(j) =    Rotational    energies.  ] ',/
     #/2x,'v  j    E(v,j; a.u.)       E(v-j_coup)     E(v; a.u.) ',
     #'    E(j; a.u.) ',/)
 670  format(///17x,'Contributions to E_vj are (Part II) : ',
     #//19x,'[  E(v,j) = E(v) + Ej1 + Ej2  ] ',/
     #/2x,'v  j    E(v,j; a.u.)        E(v;a.u.)       Ej1(a.u.) ',
     #'    Ej2(a.u.) ',/)
 680  format(///12x,'Contributions to E_vj are (Part III) : ',
     #//16x,"[   E'(v,j) = E(v) + Ej1  ] ",/
     #/2x,"v  j    E(v,j; a.u.)      E'(v,j; a.u.)     Ej1(a.u.) ",/)
 710  format(///19x,'===  Rotational energies :  === ',
     #//12x,"[   EJ(J) =  Be*[ J(J+1)-lda*lda ] + ",
     # /12x,"            D_e*[ J(J+1)-lda*lda ]**2 ",
     # /14x,"    where Be & D_e are calc. constants        ",
     #//14x,"  Er(J) = Brigid * J(J+1)                     ",
     # /14x,"     JJ = dsqrt( J(J+1) )                     ",
     # /14x,"   w(J) = dsqrt( J(J+1) )/Inertia             ",
     #//14x,"  EJ(J) = Non-rigid rotor  rot.  energies ;   ",
     # /14x,"  Er(J) = Rigid rotor rotational energies ;   ",
     # /14x,"     JJ = Classicle rot. angular momentum ;   ",
     # /14x,"   w(J) = Rotational angular velocity.      ] ",
     #//14x,'  Inertia =   u*Re*Re   = ',1PE16.8,
     # /14x,'   Brigid = 1/2*u*Re*Re = ',1PE16.8,
     #//5x,'J      EJ(J; a.u)        Er(J; a.u)         JJ  ',
     #'         w(J)   ',/)
 720  format(///9x,'[ ERRORj% = 100.0*| EJ(j) - Er(j) |/EJ(j) ] ',
     #//5x,'j      EJ(j; a.u)        Er(j; a.u)       ERRORj% ',/)
 730  format(/4x,'Average error for ROTational energies,',
     #' ERRORj% =',1PE11.4,'%',/)
 740  format(//8x,'[ ERRORj% = 100.0*| EJ(j) - Eri(j) |/EJ(j) ',
     #//8x,"   Eri(J) = B_e*[ J(J+1) - lda*lda ] + ",
     # /8x,"            D_e*[ J(J+1) - lda*lda ]**2  ",
     # /7x,"      where B_e & D_e are the INPUT constants ] ",
     #//5x,'j      EJ(j; a.u)       Eri(j; a.u)       ERRORj% ',/)
c----------------------------------------------------------
 900    return
      END
C
      subroutine  EvmaxDe4th(delc,deli,yc,yi,ev4c,ev4i)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectr0/ w0,we0,wei,wex,wey,wez,wet,wes,wer,v0,r0
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      dimension   yrc(3)
      complex*16  ycp(3)
      complex*16  y1c, y2c, y1i, y2i, z1, zi
      data  z1/(1.0, 0.0)/, zi/(0.0, 1.0)/
c----------------------------------------------------------
c    Calculate :
c  De = Evmax = We*Ymin - WeXe*Ymin**2 + WeYe*Ymin**3
c                 + WeZe*Ymin**4 ;                     (1)
c             Ymin = v_max + 1/2 .
c----------------------------------------------------------
        write( 6,100) 
        write(35,100) 
          i2 = 2
        if ( WeZe .eq. 0.0 ) i2 = 1 
      do i=1,i2
        if (i .eq. 1) then
          wx = wex
          wy = wey
          wz = wez
	  else
          wx = WeXe
          wy = WeYe
          wz = WeZe
	  endif
c----------------------------------------------------------
c   The derivative of Eq.(1) is a CUBIC equation :
c       y**3 + p*y**2 + q*y + r = 0.0
c----------------------------------------------------------
          p = 3.0*wy/(4.0*wz)
          q = - wx/(2.0*wz)
          r =   We/(4.0*wz)
C-
        call  solcubic(35,2,i,md,p,q,r,ymini,yrc,ycp)
C-
          if (i .eq. 1) yc = ymini
          if (i .eq. 2) yi = ymini
c---
      enddo
c----------------------------------------------------------
c  Find adjustment quantities delc & deli ;
c  Find the Lower & Upper bounds of solutions Ymin :
c           (y1c, y2c);     (y1i, y2i)
c
c  If deli < 0, sqrt(deli)=sqrt[(-1)*(-deli)]
c                         =sqrt(-1)*sqrt(-deli)
c                         =i*sqrt(-deli)
c----------------------------------------------------------
        delc = 9.0*wey**2  + 24.0*wez*wex
        deli = 9.0*WeYe**2 + 24.0*WeZe*WeXe
      if (delc .ge. 0.0) then
        y1c = z1*( -3.0*wey - dsqrt(delc) )/(12.0*wez)
        y2c = z1*( -3.0*wey + dsqrt(delc) )/(12.0*wez)
      else
        y1c = ( -3.0*wey - zi*dsqrt(-delc) )/(12.0*wez)
        y2c = ( -3.0*wey + zi*dsqrt(-delc) )/(12.0*wez)
      endif
        if (wez .gt. 0.0 .and. wey .lt. 0.0) then
          y1c = - y1c
          y2c = - y2c
        endif
c---
      if (deli .ge. 0.0) then
        y1i = z1*( -3.0*WeYe - dsqrt(deli) )/(12.0*WeZe)
        y2i = z1*( -3.0*WeYe + dsqrt(deli) )/(12.0*WeZe)
      else
        y1i = ( -3.0*WeYe - zi*dsqrt(-deli) )/(12.0*WeZe)
        y2i = ( -3.0*WeYe + zi*dsqrt(-deli) )/(12.0*WeZe)
      endif
        if (WeZe .gt. 0.0 .and. WeYe .lt. 0.0) then
          y1i = - y1i
          y2i = - y2i
        endif
c----------------------------------------------------------
c   Calculate the derivative of Eq.(1) and v_max.
c Since 14.43 - 0.5 = 13.93 ~=~ 14, but
c    m4c = int(14.43 - 0.50) = 13,  so we set
c    m4c = m4c + 1 
c----------------------------------------------------------
       d4c = We - 2.0*wex*yc  + 3.0*wey*yc**2  + 4.0*wez*yc**3
       d4i = We - 2.0*WeXe*yi + 3.0*WeYe*yi**2 + 4.0*WeZe*yi**3
       m4c = int(yc - 0.50)
       m4c = m4c + 1
       m4i = int(yi - 0.50)
       m4i = m4i + 1
         if ( WeZe .eq. 0.0 ) then
           d4i = 0.0
	      yi = 0.0
           m4i = 0
         endif
c----------------------------------------------------------
c  Calculate the approximate De in Eq.(1) :
c----------------------------------------------------------
      ev4c = We*yc - wex*yc**2  + wey*yc**3  + wez*yc**4
      ev4i = We*yi - WeXe*yi**2 + WeYe*yi**3 + WeZe*yi**4
c
           ky=6
         do k=1,2
           write(ky,150) 
           write(ky,160) d4c, m4c, z1*abs(y2c), yc, z1*abs(y1c)
           write(ky,170) d4i, m4i, z1*abs(y2i), yi, z1*abs(y1i)
             ky=35
         enddo 
c----------------------------------------------------------
  100 format(///12x,' === ROOTS for QUADRUPOLE form === ',
     #//6x,'Derivative of QUAdrupole form is a CUBIC equation ',
     #/9x,'            and has THREE roots : ',
     #/9x,'      [ Solved using standard method ] ',/ )
c    #/4x,'a0 > 0, there are ONE real, TWO complex conjugate ROOTS',
c    #/4x,'a0 = 0, there are THREE real ROOTS, at least TWO equal',
c    #/4x,'a0 < 0, there are THREE NON-equal real ROOTS ',/)
  150 format(//16x,'CAUTION !  You should have : ',
     #//10x,'Y2  <  Ymin  <  Y1    for  WeZe(or wez) > 0 ',
     #//10x,'Y2  >  Ymin  >  Y1    for  WeZe(or wez) < 0 ',)
  160 format(//10x,'Generated from calculated vib. constants :',
     #//11x,' d/dv {E(v_max; 4th; c)} = ',1PE12.4,
     # /9x,'Maximum vib. quantum number,  v_max(c) = ',i4,
     #//12x,' Y2_c',15x,'Ymin=Yc',16x,'Y1_c',
     #/2x,25(1H-),2x,13(1H-),2x,25(1H-),
     #/6x,"Real",7x,"Imaginary",7x,"Real",11x,"Real",7x,"Imaginary",
     #/1x,2(1PE13.5),1x,1PE14.6,1x,2(1PE13.5) )
  170 format(//10x,'Generated from the INPUT vib. constants :',
     #//11x,' d/dv {E(v_max; 4th; i)} =',1PE12.4,
     # /9x,'Maximum vib. quantum number,  v_max(i) = ',i4,
     #//12x,' Y2_i',15x,'Ymin=Yi',16x,'Y1_i',
     #/2x,25(1H-),2x,13(1H-),2x,25(1H-),
     #/6x,"Real",7x,"Imaginary",7x,"Real",11x,"Real",7x,"Imaginary",
     #/1x,2(1PE13.5),1x,1PE14.6,1x,2(1PE13.5) )
c----------------------------------------------------------
        return
      end
C
      subroutine  EvmaxDe5th(yc,yi,ev5c,ev5i)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectr0/ w0,we0,wei,wex,wey,wez,wet,wes,wer,v0,r0
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      dimension   yrc(3), xrl(15)
      complex*16  ycp(3), xrt(30)
      complex*16  ba, bm, by, ye, ea, em, x1, x2, y, z1, zi
      data  z1/(1.0, 0.0)/, zi/(0.0, 1.0)/
c----------------------------------------------------------
c    Calculate De from the QUINTIC equation :
c  De = Evmax = We*Ymin - WeXe*Ymin**2 + WeYe*Ymin**3
c                 + WeZe*Ymin**4 + WeTe*Ymin**5 ;      (1)
c             Ymin = v_max + 1/2 .
c
c  Find solutions of the equation generated from d/dv[Ev]
c  using standard method.     The derivative of Eq.(1) is
c  a equation having QUAdrupole form.
c----------------------------------------------------------
        write( 6,100) 
        write(35,100) 
          i2=2
        if ( WeTe .eq. 0.0 ) i2=1 
      do i=1,i2
        if (i .eq. 1) then
          wx = wex
          wy = wey
          wz = wez
          wt = wet
	  else
          wx = WeXe
          wy = WeYe
          wz = WeZe
          wt = WeTe
	  endif
c----------------------------------------------------------
c  Define the coefficients of the QUAdrupole equation :
c    x**4 + b*x**3 + c*x**2 + d*x + e = 0.0
c----------------------------------------------------------
           b =   4.0*wz/(5.0*wt) 
           c =   3.0*wy/(5.0*wt) 
           d = - 2.0*wx/(5.0*wt) 
           e =   We/(5.0*wt) 
c----------------------------------------------------------
c  First SOLVE the related CUBIC equation :
c        y**3 + f*y**2 + g*y + h = 0.0
c  Define the coefficients :
c----------------------------------------------------------
           f = - c
           g = b*d - 4*e
           h = - b*b*e + 4.0*c*e - d*d
c-
        call  solcubic(35,2,i,md,f,g,h,ymini,yrc,ycp)
c----------------------------------------------------------
c  Solve TWO QUAdratic equations for each REAL root y : 
c    x*x + p*x + q = 0.0 
c      p = ba or bm ;      q = ea or em .
c----------------------------------------------------------
            mk=3
          if (md .gt. 0) mk=md
        do k=1,mk
            k1 = k - 1
             y = ycp(k)
          if (md .gt. 0) y = z1*yrc(k)
            by = b*b - 4.0*c + 4.0*y
            ye = y*y - 4.0*e
          if (abs(by) .ge. 0.0) then
            ba = 0.50*z1*( b + sqrt(by) ) 
            bm = 0.50*z1*( b - sqrt(by) ) 
          else
            ba = 0.50*( b + zi*sqrt(-by) )
            bm = 0.50*( b - zi*sqrt(-by) )
          endif
c-
          if (abs(ye) .ge. 0.0) then
            ea = 0.50*z1*( y + sqrt(ye) ) 
            em = 0.50*z1*( y - sqrt(ye) ) 
          else
            ea = 0.50*( y + zi*sqrt(-ye) )
            em = 0.50*( y - zi*sqrt(-ye) )
          endif
c-
            call  quadratic(z1,ba,ea,x1,x2)
              xrt(4*k1+1) = x1
              xrt(4*k1+2) = x2
            call  quadratic(z1,bm,em,x1,x2)
              xrt(4*k1+3) = x1
              xrt(4*k1+4) = x2
        enddo
	       kk = 4*k1+4
c----------------------------------------------------------
c Find the REAL positive roots for the QUADRUPOLE equation
c----------------------------------------------------------
        if (i .eq. 1) write(35,120) 
        if (i .eq. 2) write(35,130) 
          m = 0
        do k=1,kk
            grk = real( xrt(k) )
c         if ( imag( xrt(k) ) .eq. 0.0 .and. grk .gt. 0.0 ) then 
          if ( grk .gt. 0.0 ) then 
                 m = m + 1
            xrl(m) = grk
          endif
            write(35,200) k, xrt(k)
        enddo
c---
        if ( m .gt. 0 ) then 
            if (i .eq. 1) write(35,140) 
            if (i .eq. 2) write(35,150) 
          do k=1,m
            write(35,200) k, z1*xrl(k)
          enddo
        endif
c----------------------------------------------------------
c Find the x_minimum of REAL positive roots x's.
c----------------------------------------------------------
          xmini = xrl(1)
        do k=1,m
          if ( k .gt. 1 .and. xrl(k) .lt. xrl(k-1) ) xmini = xrl(k)
        enddo
          if (i .eq. 1) yc = xmini
          if (i .eq. 2) yi = xmini
c---
	enddo
c----------------------------------------------------------
c  Calculate the derivative of Eq.(1) and v_max
c    m5c = m5c + 1 ==> Add one to get CORRECT results.
c----------------------------------------------------------
       d5c = We - 2.0*wex*yc  + 3.0*wey*yc**2  + 4.0*wez*yc**3
       d5c = d5c + 5.0*wet*yc**4 
       d5i = We - 2.0*WeXe*yi + 3.0*WeYe*yi**2 + 4.0*WeZe*yi**3
       d5i = d5i + 5.0*WeTe*yi**4 
       m5c = int(yc - 0.50)
       m5c = m5c + 1
       m5i = int(yi - 0.50)
       m5i = m5i + 1
         if ( WeTe .eq. 0.0 ) then
           d5i = 0.0
	      yi = 0.0
           m5i = 0
         endif
c==========================================================
c  Calculate the approximate De in Eq.(1) :
c----------------------------------------------------------
      ev5c=We*yc-wex*yc**2 +wey*yc**3 +wez*yc**4 +wet*yc**5
      ev5i=We*yi-WeXe*yi**2+WeYe*yi**3+WeZe*yi**4+WeTe*yi**5
        ky=6
      do k=1,2
                       write(ky,230) d5c, m5c
        if (i2 .gt. 1) write(ky,240) d5i, m5i
                       write(ky,250)
        ky=35
      enddo
c----------------------------------------------------------
  100 format(///14x,' === ROOTS for QUINTIC form === ',
     #//6x,'Derivative of QUINTIC form is a QUAdrupole equation ',
     #/9x,'            and has FOUR roots : ',
     #/9x,'      [ Solved using standard method ] ',/
     #/9x,'*** The MIDDLE results from CUBIC equation *** ',)
  120 format(/9x,'*** ',38(1H-),' ***',
     #//10x,'  Generated from calculated vib. constants :',
     #//15x,'The ROOTS of derivative equation of ',
     # /15x,'the QUINTIC vibrational energy form : ',
     # /13x,'( If they are complex, take the one with ',
     # /13x,'  minimum positive REAL part as x_mini.  ',
     #//13x,'  Each root of above CUBIC form produces ',
     # /13x,'  FOUR roots for the QUAdrupole equation. )',
     #//8x,'        m              Root_m(c)    ',
     # /15x,3(1H-),2x,31(1H-),)
  130 format(//14x,'Generated from INPUT vib. constants :',
     #//15x,'The ROOTS of derivative equation of ',
     # /15x,'the QUINTIC vibrational energy form : ',
     # /13x,'( If they are complex, take the one with ',
     # /13x,'  minimum positive REAL part as x_mini.  ',
     #//13x,'  Each root of above CUBIC form produces ',
     # /13x,'  FOUR roots for the QUAdrupole equation. )',
     #//8x,'        m              Root_m(i)    ',
     # /15x,3(1H-),2x,31(1H-),)
  140 format(//15x,'The REAL positive parts of the ROOTS ',
     #//8x,'        m              Root_m(c)    ',
     # /15x,3(1H-),2x,31(1H-),)
  150 format(//15x,'The REAL positive parts of the ROOTS ',
     #//8x,'        m              Root_m(i)    ',
     # /15x,3(1H-),2x,31(1H-),)
  200 format(15x,i2,2x,2(1PE16.8),)
  230 format(//12x,' d/dv {E(v_max; 5th; c)} = ',1PE12.4,
     #/13x,'Max. vib. quantum number, v_max(c) =',i3,)
  240 format(//12x,' d/dv {E(v_max; 5th; i)} = ',1PE12.4,
     #/13x,'Max. vib. quantum number, v_max(i) =',i3,)
  250 format(/9x,'===',7x,'Finished  for QUINTIC form',7x,'===',)
c----------------------------------------------------------
        return
      end
C
      subroutine  solcubic(kp,n,i,m,p,q,r,ymin,yrc,ycp)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      dimension   yrc(3)
      complex*16  ab, a1, b1, a2, b2
      complex*16  x1, x2, x3, z1, zi, ycp(3)
      data  z1/(1.0, 0.0)/, zi/(0.0, 1.0)/
c----------------------------------------------------------
c  SOLVE the CUBIC equation :
c        y**3 + p*y**2 + q*y + r = 0.0
c
c kp --- The channel number for printing file "fort.kp".  
c
c n = 1, Solve and print all solutions;
c   = 2, Find the positive REAL and minimum solution.
c
c i = 1, For the calculated vibrational constants;
c   = 2, For the  INPUTTED  vibrational constants.
c
c m ---  The size of the REAL solution array yrc.
c
c (p, q, r) are the coefficients of the CUBIC equation.
c   ymin -- the REAL and minimum solution.
c yrc(3) -- A  REAL   array containing REAL solutions.
c ycp(3) -- A complex array containing ALL  solutions.
c----------------------------------------------------------
c
c     As a0 < 0.0,  dsqrt(a0)=dsqrt[(-1)*(-a0)]
c                            =dsqrt(-1)*dsqrt(-a0)
c                            = i*dsqrt(-a0)
c                   dsqrt(-3)= i*dsqrt(3)
c                   dsqrt(3) = 1.732050808
c----------------------------------------------------------
       do k=1,3
         yrc(k) =  0.0
         ycp(k) = (0.0, 0.0)
	 enddo
c----------------------------------------------------------
       s3 = 1.732050808d0
        a = (3.0*q - p*p)/3.0
        b = (2.0*p*p*p - 9.0*p*q + 27.0*r)/27.0
       a0 = b*b/4.0 + a*a*a/27.0
c
      if (a0 .ge. 0.0) then
        ab = z1 * dsqrt( a0 )
      else
        ab = zi * dsqrt( - a0 )
      endif
c
c       a1 = ( - z1*b/2.0 + ab )**(1.0/3.0)
c       b1 = ( - z1*b/2.0 - ab )**(1.0/3.0)
c----------------------------------------------------------
c (-d)**1/3 = [(-1)*d]**1/3 = (-1)**1/3 * d**1/3
c           = [(-1)**1/2]**2/3 * d**1/3 
c           = i**2/3 * d**1/3
c----------------------------------------------------------
        a2 = - z1*b/2.0 + ab 
      if (imag(a2) .eq. 0.0) then
        if (real(a2) .ge. 0.0) then
          a1 = (a2)**(1.0/3.0)
        else
          a3 = abs(a2)
          a1 = (zi**(2.0/3.0) ) * (a3)**(1.0/3.0)
        endif
      else
          a1 = (a2)**(1.0/3.0)
      endif
c-
        b2 = - z1*b/2.0 - ab 
      if (imag(b2) .eq. 0.0) then
        if (real(b2) .ge. 0.0) then
          b1 = (b2)**(1.0/3.0)
        else
          b3 = abs(b2)
          b1 = (zi**(2.0/3.0) ) * (b3)**(1.0/3.0)
        endif
      else
          b1 = (b2)**(1.0/3.0)
      endif
c-
        x1 = a1 + b1
        x2 = - (a1 + b1)/2.0 + s3*zi*(a1 - b1)/2.0
        x3 = - (a1 + b1)/2.0 - s3*zi*(a1 - b1)/2.0
        ycp(1) = x1 - p/3.0
        ycp(2) = x2 - p/3.0
        ycp(3) = x3 - p/3.0
c----------------------------------------
c  Find the REAL minimum solution Ymin :
c----------------------------------------
C         write( 6,100)
C         write(kp,100)
            m = 0
c         if (i .eq. 1) write( 6,110) z1*a0,-z1*b/2.0,ab,a2,b2,a1,b1
c         if (i .eq. 1) write(kp,110) z1*a0,-z1*b/2.0,ab,a2,b2,a1,b1
c         if (i .eq. 2) write(kp,120) z1*a0,-z1*b/2.0,ab,a2,b2,a1,b1
c         if (i .eq. 2) write( 6,120) z1*a0,-z1*b/2.0,ab,a2,b2,a1,b1
C-
          if (i .eq. 1) write( 6,110)
          if (i .eq. 1) write(kp,110)
          if (i .eq. 2) write(kp,120)
          if (i .eq. 2) write( 6,120) 
        do k=1,3
            grk = real( ycp(k) )
            gik = imag( ycp(k) )
          if ( gik .eq. 0.0 ) then
                 m = m + 1
            yrc(m) = grk
          endif
            write( 6,130) k, ycp(k)
            write(kp,130) k, ycp(k)
        enddo
c-
          IF (n .eq. 1) GOTO 900
c-
        if (m .gt. 0) then
            ymin = yrc(1)
          do k=1,m
            if ( k .gt. 1 .and. yrc(k) .lt. yrc(k-1) ) ymin = yrc(k)
          enddo
        endif
c----------------------------------------------------------
  100 format(///10x,' === ROOT range for QUADRUPOLE form === ',
     #//6x,'Derivative of QUAdrupole form is a CUBIC equation ',
     #/9x,'            and has THREE roots : ',/)
c    #/4x,'a0 > 0, there are ONE real, TWO complex conjugate ROOTS',
c    #/4x,'a0 = 0, there are THREE real ROOTS, at least TWO equal',
c    #/4x,'a0 < 0, there are THREE NON-equal real ROOTS ',/)
  110 format(/14x,'From calculated vib. constants :',
C 110 format(/10x,'CHECKS from calculated vib. constants :',
c    #/26x,'ab = sqrt(a0) ',
c    #/15x,'a0 = ',2(1PE12.4),
c    #/13x,'-b/2 = ',2(1PE12.4),
c    #/15x,'ab = ',2(1PE12.4),/
c    #/13x,'a2 = -b/2 + ab ;   b2 = -b/2 - ab ',
c    #/15x,'a2 = ',2(1PE12.4),
c    #/15x,'b2 = ',2(1PE12.4),/
c    #/13x,'a1 = (a2)**(1/3) ; b1 = (b2)**(1/3) ',
c    #/15x,'a1 = ',2(1PE12.4),
c    #/15x,'b1 = ',2(1PE12.4),/
     #//7x,'        Y_mini = Y_root = v_max + 1/2  ',
     # /7x,'Y_mini = The REAL positive root with minimum value ',
     #/10x,'      n            Root_n(c)  ',
     #/15x,3(1H-),3x,25(1H-),)
  120 format(//14x,'From the INPUT  vib. constants :',
C 120 format(//10x,'CHECK from the INPUT  vib. constants :',
c    #/26x,'ab = sqrt(a0) ',
c    #/15x,'a0 = ',2(1PE12.4),
c    #/13x,'-b/2 = ',2(1PE12.4),
c    #/15x,'ab = ',2(1PE12.4),/
c    #/13x,'a2 = -b/2 + ab ;   b2 = -b/2 - ab ',
c    #/15x,'a2 = ',2(1PE12.4),
c    #/15x,'b2 = ',2(1PE12.4),/
c    #/13x,'a1 = (a2)**(1/3) ; b1 = (b2)**(1/3) ',
c    #/15x,'a1 = ',2(1PE12.4),
c    #/15x,'b1 = ',2(1PE12.4),/
     #//7x,'Y_mini = The REAL positive root with minimum value ',
     #/10x,'      n            Root_n(i)  ',
     #/15x,3(1H-),3x,25(1H-),)
  130 format(15x,i2,3x,2(1PE13.5) )
c----------------------------------------------------------
  900   return
      end
C
      subroutine quadratic(a,b,c,x1,x2)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      complex*16  a, b, c, p, q, x1, x2
c----------------------------------------------------------
C   Find the roots, x1 & x2, of a quadratic equation
C               a*x*x + b*x + c = 0
C using the skills on P. 178 of "Numerical Reccipes".
c       write(6,*) '  sign(1.0, b) = ', sign(1.0,b)
c----------------------------------------------------------
c     q = -0.5d0*(b + sign(1.0,b)*dsqrt( b*b - 4.0*a*c) )
c       x1 = q/a
c       x2 = c/q
c-
      p = 0.5d0*(- b + sqrt( b*b - 4.0*a*c ) )
      q = 0.5d0*(- b - sqrt( b*b - 4.0*a*c ) )
        x1 = p/a
        x2 = q/a
c---
        return
      end
C
C
      subroutine  vjcoef(We1)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectrb/ Be,alphae,gamae,Der,betae,bryd,wvib0(40)
      common /vjcoef0/ a0,a1,a2,a3,a4,a5,a6,a7,a8 
      common /vjcoef1/ b1,b3,b4,b5,b6,b7,b8
c----------------------------------------------------------
c   Calculate the coefficients used to determine the
c vibrational-rotational constants WeXe,..., Alfa_e,...
c----------------------------------------------------------
       a = amu*We1
      ar = dsqrt(a)
      a0 =   1.0/(Re*Re)   
      a1 = - 2.0/(Re*Re*Re*ar)   
      a2 =   3.0/(Re*Re*Re*Re*a)   
      a3 = - 4.0/( Re**5 * ar**3)   
      a4 =   5.0/( Re**6 *  a**2)   
      a5 = - 6.0/( Re**7 * ar**5)   
      a6 =   7.0/( Re**8 *  a**3)   
      a7 = - 8.0/( Re**9 * ar**7)   
      a8 =   9.0/(Re**10 *  a**4)
c
      b1 =   1.0/(             ar)
      b3 =   1.0/(    6.0 * ar**3)
      b4 =   1.0/(   24.0 *  a**2)
      b5 =   1.0/(  120.0 * ar**5)
      b6 =   1.0/(  720.0 *  a**3)
      b7 =   1.0/( 5040.0 * ar**7)
      b8 =   1.0/(40320.0 *  a**4)
c----------------------------------------------------------
        return
      END
C
      subroutine  convj(gg,nd,We1)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectrb/ Be,alphae,gamae,Der,betae,bryd,wvib0(40)
      common /spectr0/ w0,we0,wei,wex,wey,wez,wet,wes,wer,v0,r0
      common /spectr1/ Bee,ale,gae,eta3,eta4,eta5,eta6,eta7
      common /spectr2/ Dee,bete,xsi2,xsi3,xsi4,xsi5,xsi6,xsi7
      common /vjcoef0/ a0,a1,a2,a3,a4,a5,a6,a7,a8 
      common /vjcoef1/ b1,b3,b4,b5,b6,b7,b8
      common /forcecon/ f1,f2,f3,f4,f5,f6,f7,f8
      common /forceco1/ f9,f10,f11,f12,f13,f14,f15
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      dimension  gg(nd)
c----------------------------------------------------------
      f1 = gg(1)     
      f2 = gg(2)     
      f3 = gg(3)     
      f4 = gg(4)     
      f5 = gg(5)     
      f6 = gg(6)     
      f7 = gg(7)     
      f8 = gg(8)     
c---
      f12 = f1*f1
      f32 = f3*f3
      f42 = f4*f4
      f52 = f5*f5
      f62 = f6*f6
      f72 = f7*f7
      f82 = f8*f8
c
      a12 = a1*a1
      a22 = a2*a2
      a32 = a3*a3
      a42 = a4*a4
      a52 = a5*a5
      a62 = a6*a6
      a72 = a7*a7
      a82 = a8*a8
c
      b12 = b1*b1
      b32 = b3*b3
      b42 = b4*b4
      b52 = b5*b5
      b62 = b6*b6
      b72 = b7*b7
      b82 = b8*b8
c
      w0 = 3.0*b4*f4/8.0 + 315.0*b8*f8/128.0 - 7.0*b32*f32/(16.0*We1)
      w0 = w0 - 1107.0*b52*f52/(256.0*We1) -945.*b4*b6*f4*f6/(128.*We1)
      w0 = w0 - 1155.0*b3*b7*f3*f7/(128.0*We1)
      w0 = w0 - 180675.0*b72*f72/(2048.0*We1)
      w0 = w0 - 89775.0*b6*b8*f6*f8/(512.0*We1)
c
c-------------------------------------------------------
c  If We1 = We + we0, then we0 is a function of itself.
c-------------------------------------------------------
      we0 = 25.0*b6*f6/8.0 - 67.0*b42*f42/(16.0*We1)
      we0 = we0 - 95.0*b3*b5*f3*f5/(8.0*We1) -19277.*b62*f62/(256.*We1)
      we0 = we0 - 22029.0*b5*b7*f5*f7/(128.0*We1)
      we0 = we0 - 10521.0*b4*b8*f4*f8/(64.0*We1)
      we0 = we0 - 5450499.0*b82*f82/(2048.0*We1)
c
      wex = 3.0*b4*f4/2.0 + 245.0*b8*f8/16.0 - 15.0*b32*f32/(4.0*We1)
      wex = wex - 1085.*b52*f52/(32.0*We1) -885.*b4*b6*f4*f6/(16.*We1)
      wex = wex - 1365.0*b3*b7*f3*f7/(16.0*We1)
      wex = wex - 444381.0*b72*f72/(512.0*We1)
      wex = wex - 204771.0*b6*b8*f6*f8/(128.0*We1)
        if (meny .eq. 0) wex = - wex 
c
      wey = 5.0*b6*f6/2.0 - 17.0*b42*f42/(4.0*We1) 
      wey = wey - 35.0*b3*b5*f3*f5/(2.*We1) -4145.0*b62*f62/(32.0*We1)
      wey = wey - 5355.0*b5*b7*f5*f7/(16.0*We1)
      wey = wey - 2205.0*b4*b8*f4*f8/(8.0*We1)
      wey = wey - 2947595.0*b82*f82/(512.0*We1)
c
      wez = 35.0*b8*f8/8.0 - 315.*b52*f52/(16.*We1) 
      wez = wez - 165.*b4*b6*f4*f6/(8.*We1) -315.*b3*b7*f3*f7/(8.*We1)
      wez = wez - 82005.*b72*f72/(128.*We1) 
      wez = wez - 33845.*b6*b8*f6*f8/(32.*We1)
C       if (meny .eq. 0) wez = - wez 
c
      wet = - 393.0*b62*f62/(16.0*We1) - 693.0*b5*b7*f5*f7/(8.0*We1)
      wet = wet - 189.*b4*b8*f4*f8/(4.*We1) 
      wet = wet - 239841.*b82*f82/(128.0*We1)
c
      wes = - 3003.0*b72*f72/(32.0*We1) - 889.0*b6*b8*f6*f8/(8.0*We1)
C       if (meny .eq. 0) wes = - wes 
c
      wer = - 3985.0*b82*f82/(32.0*We1)
c---
      Bee = a0 + 3.0*a4/8.0 + 315.0*a8/128.0 - 7.0*a3*b3*f3/(8.0*We1)     
      Bee = Bee - 1155.0*a7*b3*f3/(128.0*We1) -3.0*a2*b4*f4/(4.0*We1)
      Bee = Bee - 945.0*a6*b4*f4/(128.0*We1) -15.0*a1*b5*f5/(8.0*We1)
      Bee = Bee - 1107.*a5*b5*f5/(128.0*We1) -945.*a4*b6*f6/(128.*We1)
      Bee = Bee - 89775.*a8*b6*f6/(512.*We1)-1155.*a3*b7*f7/(128.*We1)
      Bee = Bee - 180675.*a7*b7*f7/(1024.*We1)-315.*a2*b8*f8/(32.*We1)
      Bee = Bee - 89775.0*a6*b8*f8/(512.0*We1)
      Bee = Bee/(2.0*amu)
c
      ale = a2 + 25.0*a6/8.0 - 3.0*a1*b3*f3/We1 -95.0*a5*b3*f3/(8.*We1)
      ale = ale - 67.0*a4*b4*f4/(8.0*We1) - 10521.0*a8*b4*f4/(64.0*We1)
      ale = ale - 95.0*a3*b5*f5/(8.0*We1) - 22029.0*a7*b5*f5/(128.*We1)
      ale = ale - 75.0*a2*b6*f6/(8.0*We1) - 19277.0*a6*b6*f6/(128.*We1)
      ale = ale - 175.0*a1*b7*f7/(8.*We1) - 22029.0*a5*b7*f7/(128.*We1)
      ale = ale - 10521.0*a4*b8*f8/(64.0*We1) 
      ale = ale - 5450499.0*a8*b8*f8/(1024.0*We1)
      ale = - ale/(2.0*amu)
c
      gae = 3.0*a4/2.0 + 245.0*a8/16.0 - 15.0*a3*b3*f3/(2.0*We1)
      gae = gae - 1365.0*a7*b3*f3/(16.0*We1) - 3.0*a2*b4*f4/We1
      gae = gae - 885.0*a6*b4*f4/(16.0*We1) - 15.0*a1*b5*f5/(2.0*We1)
      gae = gae - 1085.0*a5*b5*f5/(16.*We1) - 885.*a4*b6*f6/(16.0*We1)
      gae = gae - 204771.0*a8*b6*f6/(128.0*We1) 
      gae = gae - 1365.*a3*b7*f7/(16.*We1) -444381.*a7*b7*f7/(256.*We1) 
      gae = gae - 245.0*a2*b8*f8/(4.0*We1) -204771.*a6*b8*f8/(128.*We1) 
      gae = gae/(2.0*amu)
c
c--- Bee, ale, gae   are calculated by DEFINITION !
c
      et3 = 5.0*a6/2.0 - 35.0*a5*b3*f3/(2.0*We1) -17.*a4*b4*f4/(2.*We1)
      et3 = et3 - 2205.0*a8*b4*f4/(8.0*We1) - 35.0*a3*b5*f5/(2.0*We1)
      et3 = et3 - 5355.0*a7*b5*f5/(16.0*We1) - 15.*a2*b6*f6/(2.0*We1)
      et3 = et3 - 4145.0*a6*b6*f6/(16.0*We1) - 35.*a1*b7*f7/(2.0*We1)
      et3 = et3 - 5355.0*a5*b7*f7/(16.0*We1) - 2205.*a4*b8*f8/(8.*We1)
      et3 = et3 - 2947595.0*a8*b8*f8/(256.0*We1)
      eta3 = et3/(2.0*amu)
c
      et4 = 35.*a8/8.0 - 315.*a7*b3*f3/(8.*We1) -165.*a6*b4*f4/(8.*We1)
      et4 = et4 - 315.0*a5*b5*f5/(8.0*We1) - 165.0*a4*b6*f6/(8.0*We1)
      et4 = et4 - 33845.0*a8*b6*f6/(32.0*We1) -315.0*a3*b7*f7/(8.0*We1)
      et4 = et4 - 82005.0*a7*b7*f7/(64.0*We1) - 35.0*a2*b8*f8/(2.0*We1)
      et4 = et4 - 33845.0*a6*b8*f8/(32.0*We1)
      eta4 = et4/(2.0*amu)
c
      et5 = - 189.0*a8*b4*f4/(4.0*We1) - 693.0*a7*b5*f5/(8.0*We1)
      et5 = et5 - 393.0*a6*b6*f6/(8.0*We1) - 693.0*a5*b7*f7/(8.0*We1)
      et5 = et5 - 189.0*a4*b8*f8/(4.0*We1) - 239841.*a8*b8*f8/(64.*We1)
      eta5 = et5/(2.0*amu)
c
      et6 = - 889.0*a8*b6*f6/(8.0*We1) - 3003.0*a7*b7*f7/(16.0*We1)
      et6 = et6 - 889.0*a6*b8*f8/(8.0*We1)
      eta6 = et6/(2.0*amu)
c
      eta7 = - 3985.0*a8*b8*f8/(32.0*amu*We1)
c---
      Dee = a12/(2.0*We1) + 7.0*a32/(16.0*We1) + 3.0*a2*a4/(4.0*We1)
      Dee = Dee + 15.0*a1*a5/(8.0*We1) + 1107.0*a52/(256.0*We1)
      Dee = Dee + 945.0*a4*a6/(128.0*We1) + 1155.0*a3*a7/(128.0*We1)
      Dee = Dee + 180675.0*a72/(2048.0*We1) + 315.0*a2*a8/(32.0*We1)
      Dee = Dee + 89775.0*a6*a8/(512.0*We1)
      Dee = Dee/(4.0*amu*amu)
c
      bet = a22/(2.0*We1) + 3.0*a1*a3/We1 + 67.0*a42/(16.0*We1)
      bet = bet + 95.0*a3*a5/(8.0*We1) + 75.0*a2*a6/(8.0*We1)
      bet = bet + 19277.0*a62/(256.0*We1) + 175.0*a1*a7/(8.0*We1)
      bet = bet + 22029.0*a5*a7/(128.0*We1) + 10521.0*a4*a8/(64.0*We1)
      bet = bet + 5450499.0*a82/(2048.0*We1)
      bete = bet/(4.0*amu*amu)
c
      xs2 = 15.0*a32/(4.0*We1) + 3.0*a2*a4/We1 + 15.0*a1*a5/(2.0*We1)
      xs2 = xs2 + 1085.0*a52/(32.0*We1) + 885.0*a4*a6/(16.0*We1)
      xs2 = xs2 + 1365.0*a3*a7/(16.0*We1) + 444381.0*a72/(512.0*We1)
      xs2 = xs2 + 245.0*a2*a8/(4.0*We1) + 204771.0*a6*a8/(128.0*We1)
      xsi2 = xs2/(4.0*amu*amu)
c
      xs3 = 17.0*a42/(4.0*We1) +35.0*a3*a5/(2.*We1) +15.*a2*a6/(2.*We1)
      xs3 = xs3 + 4145.0*a62/(32.0*We1) + 35.0*a1*a7/(2.0*We1)
      xs3 = xs3 + 5355.0*a5*a7/(16.0*We1) + 2205.0*a4*a8/(8.0*We1)
      xs3 = xs3 + 2947595.0*a82/(512.0*We1)
      xsi3 = xs3/(4.0*amu*amu)
c
      xs4 = 315.0*a52/(16.0*We1) + 165.0*a4*a6/(8.0*We1)
      xs4 = xs4 + 315.0*a3*a7/(8.0*We1) + 82005.0*a72/(128.0*We1)
      xs4 = xs4 + 35.0*a2*a8/(2.0*We1) + 33845.0*a6*a8/(32.0*We1)
      xsi4 = xs4/(4.0*amu*amu)
c
      xs5 = 393.0*a62/(16.0*We1) + 693.0*a5*a7/(8.0*We1)
      xs5 = xs5 + 189.0*a4*a8/(4.0*We1) + 239841.0*a82/(128.0*We1)
      xsi5 = xs5/(4.0*amu*amu)
c
      xs6 = 3003.0*a72/(32.0*We1) + 889.0*a6*a8/(8.0*We1)
      xsi6 = xs6/(4.0*amu*amu)
c
      xsi7 = 3985.0*a82/(128.0*amu*amu*We1)
c
c----------------------------------------------------------
        return
      END
C===
      subroutine  vibconst(wef,fns,wvib,kv)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectr0/ w0,we0,wei,wex,wey,wez,wet,wes,wer,v0,r0
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      dimension   wvib(kv),fns(kv)
c----------------------------------------------------------
        a = amu*wef
       ar = dsqrt(a)
       b1 =   1.0/(             ar)
       b3 =   1.0/(    6.0 * ar**3)
       b4 =   1.0/(   24.0 *  a**2)
       b5 =   1.0/(  120.0 * ar**5)
       b6 =   1.0/(  720.0 *  a**3)
       b7 =   1.0/( 5040.0 * ar**7)
       b8 =   1.0/(40320.0 *  a**4)
c-
      b32 = b3*b3
      b42 = b4*b4
      b52 = b5*b5
      b62 = b6*b6
      b72 = b7*b7
      b82 = b8*b8
c---
      f3 = fns(3)     
      f4 = fns(4)     
      f5 = fns(5)     
      f6 = fns(6)     
      f7 = fns(7)     
      f8 = fns(8)     
c-
      f32 = f3*f3
      f42 = f4*f4
      f52 = f5*f5
      f62 = f6*f6
      f72 = f7*f7
      f82 = f8*f8
c----------------------------------------------------------
      w0 = 3.0*b4*f4/8.0 + 315.0*b8*f8/128.0 - 7.0*b32*f32/(16.0*wef)
      w0 = w0 - 1107.0*b52*f52/(256.0*wef) -945.*b4*b6*f4*f6/(128.*wef)
      w0 = w0 - 1155.0*b3*b7*f3*f7/(128.0*wef)
      w0 = w0 - 180675.0*b72*f72/(2048.0*wef)
      w0 = w0 - 89775.0*b6*b8*f6*f8/(512.0*wef)
c
Cc    w0 = - w0
c
c-------------------------------------------------------
c  If wef = We + we0, then we0 is a function of itself.
c-------------------------------------------------------
      we0 = 25.0*b6*f6/8.0 - 67.0*b42*f42/(16.0*wef)
      we0 = we0 - 95.0*b3*b5*f3*f5/(8.0*wef) -19277.*b62*f62/(256.*wef)
      we0 = we0 - 22029.0*b5*b7*f5*f7/(128.0*wef)
      we0 = we0 - 10521.0*b4*b8*f4*f8/(64.0*wef)
      we0 = we0 - 5450499.0*b82*f82/(2048.0*wef)
c
Cc    we0 = - we0 
c
        wvib(30) = w0
        wvib(31) = we0
c----------------------------------------------------------
      if (nw0 .ge. 1) wvib(2) = wef 
      if (nw0 .lt. 1) wvib(1) = wef 
c-
      wex = 3.0*b4*f4/2.0 + 245.0*b8*f8/16.0 - 15.0*b32*f32/(4.0*wef)
      wex = wex - 1085.*b52*f52/(32.0*wef) -885.*b4*b6*f4*f6/(16.*wef)
      wex = wex - 1365.0*b3*b7*f3*f7/(16.0*wef)
      wex = wex - 444381.0*b72*f72/(512.0*wef)
      wex = wex - 204771.0*b6*b8*f6*f8/(128.0*wef)
        if (meny .eq. 0)  wex = - wex 
      if (nw0 .ge. 1) wvib(3) =   wex 
      if (nw0 .lt. 1) wvib(2) =   wex 
c
      wey = 5.0*b6*f6/2.0 - 17.0*b42*f42/(4.0*wef) 
      wey = wey - 35.0*b3*b5*f3*f5/(2.*wef) -4145.0*b62*f62/(32.0*wef)
      wey = wey - 5355.0*b5*b7*f5*f7/(16.0*wef)
      wey = wey - 2205.0*b4*b8*f4*f8/(8.0*wef)
      wey = wey - 2947595.0*b82*f82/(512.0*wef)
        if (nw0 .ge. 1) wvib(4) =   wey 
        if (nw0 .lt. 1) wvib(3) =   wey 
c
      wez = 35.0*b8*f8/8.0 - 315.*b52*f52/(16.*wef) 
      wez = wez - 165.*b4*b6*f4*f6/(8.*wef) -315.*b3*b7*f3*f7/(8.*wef)
      wez = wez - 82005.*b72*f72/(128.*wef) 
      wez = wez - 33845.*b6*b8*f6*f8/(32.*wef)
C       if (meny .eq. 0)  wez = - wez 
      if (nw0 .ge. 1) wvib(5) =   wez 
      if (nw0 .lt. 1) wvib(4) =   wez 
c
      wet = - 393.0*b62*f62/(16.0*wef) - 693.0*b5*b7*f5*f7/(8.0*wef)
      wet = wet - 189.*b4*b8*f4*f8/(4.*wef) 
      wet = wet - 239841.*b82*f82/(128.0*wef)
        if (nw0 .ge. 1) wvib(6) =   wet 
        if (nw0 .lt. 1) wvib(5) =   wet 
c
      wes = - 3003.0*b72*f72/(32.0*wef) - 889.0*b6*b8*f6*f8/(8.0*wef)
C       if (meny .eq. 0)  wes = - wes 
      if (nw0 .ge. 1) wvib(7) = wes 
      if (nw0 .lt. 1) wvib(6) = wes 
c
      wer = - 3985.0*b82*f82/(32.0*wef)
        if (nw0 .ge. 1) wvib(8) =   wer 
        if (nw0 .lt. 1) wvib(7) =   wer 
c
c----------------------------------------------------------
        return
      END
c===
      subroutine  funcv(np,n,x,fvec)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      dimension   fvec(np), x(np)
      common /vjcoef1/ b1,b3,b4,b5,b6,b7,b8
      common /Eswitc4/ cxe,cye,cze,cte,cse,cre
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /fncom/ xmean(30),gf(30)
C-------------------------------------------------------------    
C  As passed in, x(n) contains INITIAL guess values for Fn's.
C-------------------------------------------------------------    
      if (nw0 .ge. 1) then
         Wei = xmean(2)
        wexi = xmean(3)
        weyi = xmean(4)
        wezi = xmean(5)
        weti = xmean(6)
        wesi = xmean(7)
        weri = xmean(8)
      else
         Wei = xmean(1)
        wexi = xmean(2)
        weyi = xmean(3)
        wezi = xmean(4)
        weti = xmean(5)
        wesi = xmean(6)
        weri = xmean(7)
	endif
C------------------------------
      if (nlin .eq. 1) wexi = wexi*cxe
c----------------------------------------------------------
        call vjcoef(Wei)
c----------------------------------------------------------
      b12 = b1*b1
      b32 = b3*b3
      b42 = b4*b4
      b52 = b5*b5
      b62 = b6*b6
      b72 = b7*b7
      b82 = b8*b8
c----------------------------------------------------------
       x1 = 3.0*b4/2.0
       x2 = 245.0*b8/16.0
       x3 = 15.0*b32/(4.0*Wei)
       x4 = 1085.*b52/(32.0*Wei)
       x5 = 885.*b4*b6/(16.*Wei)
       x6 = 1365.0*b3*b7/(16.0*Wei)
       x7 = 444381.0*b72/(512.0*Wei)
       x8 = 204771.0*b6*b8/(128.0*Wei)
c
       y1 = 5.0*b6/2.0
       y2 = 17.0*b42/(4.0*Wei)
       y3 = 35.0*b3*b5/(2.*Wei)
       y4 = 4145.0*b62/(32.0*Wei)
       y5 = 5355.0*b5*b7/(16.0*Wei)
       y6 = 2205.0*b4*b8/(8.0*Wei)
       y7 = 2947595.0*b82/(512.0*Wei)
c
       z1 = 35.0*b8/8.0
       z2 = 315.*b52/(16.*Wei)
       z3 = 165.*b4*b6/(8.*Wei)
       z4 = 315.*b3*b7/(8.*Wei)
       z5 = 82005.*b72/(128.*Wei)
       z6 = 33845.*b6*b8/(32.*Wei)
c
       t1 = 393.0*b62/(16.0*Wei)
       t2 = 693.0*b5*b7/(8.0*Wei)
       t3 = 189.*b4*b8/(4.*Wei)
       t4 = 239841.*b82/(128.0*Wei)
c
       s1 = 3003.0*b72/(32.0*Wei)
       s2 = 889.0*b6*b8/(8.0*Wei)
c----------------------------------------------------------
c  x(1)=f3; x(2)=f4; x(3)=f5; x(4)=f6; x(5)=f7; x(6)=f8.
c-
c     fvec(1) = wexi + x1*x(2) + x2*x(6) - x3*x(1)*x(1)
c    # - x4*x(3)*x(3) - x5*x(2)*x(4) - x6*x(1)*x(5)
c    # - x7*x(5)*x(5) - x8*x(4)*x(6)
c----------------------------------------------------------
      fvec(1) = - wexi + x3*x(1)*x(1) - x1*x(2) - x2*x(6)
Cc    fvec(1) =   wexi + x3*x(1)*x(1) - x1*x(2) - x2*x(6)
     # + x4*x(3)*x(3) + x5*x(2)*x(4) + x6*x(1)*x(5)
     # + x7*x(5)*x(5) + x8*x(4)*x(6)
c
      fvec(2) = - weyi + y1*x(4) - y2*x(2)*x(2)
Cc    fvec(2) =   weyi + y1*x(4) - y2*x(2)*x(2)
     # - y3*x(1)*x(3) - y4*x(4)*x(4) - y5*x(3)*x(5)
     # - y6*x(2)*x(6) - y7*x(6)*x(6) 
c
      fvec(3) = - wezi + z1*x(6) - z2*x(3)*x(3) - z3*x(2)*x(4) 
Cc    fvec(3) =   wezi + z1*x(6) - z2*x(3)*x(3) - z3*x(2)*x(4) 
     # - z4*x(1)*x(5) - z5*x(5)*x(5) - z6*x(4)*x(6) 
c
      fvec(4) = - weti - t1*x(4)*x(4) - t2*x(3)*x(5) 
Cc    fvec(4) =   weti - t1*x(4)*x(4) - t2*x(3)*x(5) 
     # - t3*x(2)*x(6) - t4*x(6)*x(6) 
c
      fvec(5) = - wesi - s1*x(5)*x(5) - s2*x(4)*x(6) 
Cc    fvec(5) =   wesi - s1*x(5)*x(5) - s2*x(4)*x(6) 
c
      fvec(6) = - weri - x(6)*x(6) * 3985.0*b82/(32.0*Wei)
Cc    fvec(6) =   weri - x(6)*x(6) * 3985.0*b82/(32.0*Wei)
c----------------------------------------------------------
        RETURN
      END
C===
      REAL*8 FUNCTION FACX(I)
      IMPLICIT REAL*8 (A-H,O-Z)
C---------------------------------------------------------------------
C     THIS IS A FACTORIAL FUNCTION  I!
C---------------------------------------------------------------------
      DIMENSION TABLE(15)
      DATA TABLE/1.0D+00,2.0D+00,6.0D+00,24.0D+00,120.0D+00,720.0D+00,
     1      5040.0D+00,40320.0D+00,362880.0D+00,
     2      36288.0D+02,399168.0D+02,4790016.0D+02,62270208.0D+02,
     3      871782912.0D+02,1307674368.0D+03/
C---------------------------------------------------------------------
      FACX=1.0
      IF (I)95,100,10
   10 IF (I-15)20,20,30
   20 FACX=TABLE(I)
      GO TO 200
C-----------------------
   30 FJ=16.0D+00
      FACX=TABLE(15)
      DO 40 J=16,I
      FACX=FACX*FJ
   40 FJ=FJ+1.0D+00
      GO TO 200
   95 FACX=0.0D+00
  100 CONTINUE
C-----------------------
  200 RETURN
      END
C
      real*8 function Pnj (n,j)
      IMPLICIT REAL*8 (A-H,O-Z)
      if ( j .gt.n ) then
          Pnj = 0.0
      return
      endif
      Pnj = facx (n)/  facx (n - j )
      return
      end               
C
      real*8 function Cnj (n,j)
      IMPLICIT REAL*8 (A-H,O-Z)
      if ( j .gt.n ) then
      Cnj = 0.0
      return
      endif
      Cnj = facx (n)/ (facx (n-j) * facx ( j ))
      return
      end               
C===
       subroutine   FindnumDe(kr,key,nv,nvp,Decal,ErrDe)
       implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
c===================================================================    
c   TRANDITIONAL expression (for meny = 0) :   
c Ev = We(v+1/2)-WeXe(v+1/2)^2+WeYe(v+1/2)^3+WeZe(v+1/2)^4
c       +WeTe(v+1/2)^5+WeSe(v+1/2)^6+WeRe(v+1/2)^7            ...(1)
c===
c   NEW  expression (for meny = 1) :   
c Ev = We(v+1/2)-WeXe(v+1/2)^2+WeYe(v+1/2)^3-WeZe(v+1/2)^4
c       +WeTe(v+1/2)^5-WeSe(v+1/2)^6+WeRe(v+1/2)^7 - ...      ...(2)
c
c dEv
c --- = We-2WeXe(v+1/2)+3WeYe(v+1/2)^2-4WeZe(v+1/2)^3  
c  dv   +5WeTe(v+1/2)^4-6WeSe(v+1/2)^5+7WeRe(v+1/2)^6 == 0.0  ...(3)  
c
c  To find the MINIMUM POSITIVE roots for equations like Eq.(2) .
c--------------------------------------------------------------------
      dimension   CON(20),DisE(20),DisEc(20),DisEi(20),DisE0(20)
      dimension   CONc(20),CONi(20),dltDec(20),dltDei(20),DltE0(20)
      dimension   Numax(20),Numaxc(20),Numaxi(20)
      dimension   der0(20),der1(20),der2(20),der3(20) 
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectr0/ w0,we0,wei,wex,wey,wez,wet,wes,wer,v0,r0
      common /spectr3/ wen,wxe,wye,wze,wte,wse,wre
      common /spectr4/ w08,w09,w10,w11,w12,w13,w14,w15,w16
      common /Eswitch/ aye,aze,ate,ase,are
      common /Eswitc4/ cxe,cye,cze,cte,cse,cre
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /LA1/ Deinp,Ner,mset,mcst,mtd,mdp,mv,nvs(20,80)
c---
        if (kr .eq. 0) then 
          CONi(1)=We 
          CONi(2)=-WeXe
          CONi(3)=aye*WeYe
          CONi(4)= aze*WeZe
          CONi(5)=ate*WeTe
          CONi(6)= ase*WeSe
          CONi(7)=are*WeRe
c-
          CONc(1)=We + we0
          CONc(2)=-wex
          CONc(3)=aye*wey
          CONc(4)= aze*wez
          CONc(5)=ate*wet
          CONc(6)= ase*wes
          CONc(7)=are*wer
        else
          CONi(1)=wen
C         CONi(2)=-cxe*wxe
          CONi(2)= cxe*wxe
          CONi(3)=cye*wye
          CONi(4)= cze*wze
          CONi(5)=cte*wte
          CONi(6)= cse*wse
          CONi(7)=cre*wre
        endif
c---
          CON(8)= w08
          CON(9)= w09
          CON(10)=w10
          CON(11)=w11
          CON(12)=w12
          CON(13)=w13
          CON(14)=w14
          CON(15)=w15
c-----------------------------------------------
c  nn = 1, for INPUT vib. constants;
c     = 2, for Calc. vib. constants.
c-----------------------------------------------
        if (kr .eq. 0) nm=2
        if (kr .gt. 0) nm=1
c-
      DO 100 nn=1,nm
             kk=0
c-
C     DO  60  i=2,7 
      DO  60  i=2,mcst
C         endv = 100.5
C         toll = 1.0E-08
          toll = 1.0E-15
          star = 0.5
          endv = imv*1.0 + 0.5
        DO m = 1,7
          if (nn .eq. 1)then
            CON(m)=CONi(m)
          elseif(nn .eq. 2)then
            CON(m)=CONc(m)
          endif
c----------------------------------------------------------------
c   For key > 0, using the VIBrational constants from  A*X = B
c to cal. the NUMerical De's and derivatives del_E(v)/del_v. 
c Since these constants alternate signs "+" & "-" interchangely,
c we take their absolute values here, and carry their signs back
c when cal. De (--> DisE) & deriv. (--> tt, OR tt1) later. 
c----------------------------------------------------------------
          if ( key .gt. 0)  CON(m)=abs( CON(m) )
        ENDDO
c---
  10    continue
c          step = (endv - star)/100.0
           step = (endv - star)/(imv*1.0)
            x=star
c
      DO 30 k=1,imv
         x = x + step
c----------------------------------------------------------------
c  tt = CON1 - 2.0*CON2*x + 3.0*CON3*x**2 + 4.0*CON4*x**3 
c        + 5.0*CON5*x**4 + 6.0*CON6*x**5 +7.0*CON7*x**6
c
c  x1 = x  + step
c tt1 = CON1 - 2.0*CON2*x1 + 3.0*CON3*x1**2 + 4.0*CON4*x1**3 
c        + 5.0*CON5*x1**4 + 6.0*CON6*x1**5 +7.0*CON7*x1**6
c
c  Find the derivatives f_I for each I : tt
c----------------------------------------------------------------
         tt=0.0
       DO j=1,I
         if (key .eq. 0) then
           tt = tt + j * CON(j)*x**(j-1)
         else
           tt = tt + (-1.0)**(j-1) * j * CON(j)*x**(j-1)
         endif
       ENDDO
          x1 = x + step 
         tt1 = 0.0
       DO j=1,I
         if (key .eq. 0) then
           tt1 = tt1 + j * CON(j)*x1**(j-1)
         else
           tt1 = tt1 + (-1.0)**(j-1) * j * CON(j)*x1**(j-1)
         endif
       ENDDO
c----------------------------------------------------------------
c  Numaxi/Numaxc contains the maximum vib. quantum number.
c    Numax = Numax + 1 ==> Add to get CORRECT results.
c----------------------------------------------------------------
       if (abs(tt) .lt. toll) then
                y=x
c        Numax(i)= int( y - 0.50d0 )    
c        Numax(i)= int( y - 0.50d0 ) + 1    
         Numax(i)= int( y  )   
          if(nn .eq. 1)then
             Numaxi(i)=Numax(i)
               der1(i)=tt
               der3(i)=tt1
          elseif(nn .eq. 2)then
             Numaxc(i)=Numax(i)
               der0(i)=tt
               der2(i)=tt1
          endif 
             kk=1
           goto 40   
       elseif(abs(tt).gt.toll .and. tt*tt1 .lt. 0.0) then
         star=x
         endv=x1
           goto 10
       elseif(i .gt. 2 .and. nn .eq. 1) then
                y=x
c        Numax(i)= int( y - 0.50d0 )    
c        Numax(i)= int( y - 0.50d0 ) + 1    
         Numax(i)= int( y )    
             Numaxi(i)=Numax(i)
               der1(i)=tt
               der3(i)=tt1
       endif
c--
 30   CONTINUE
c--
 40   CONTINUE                       
c--
      if ( kk .ne. 1 .and. kr .eq. 0 ) then
        write( 6,280) star, endv
        write(38,280) star, endv
      endif
c----------------------------------------------------------------
c   Calculate the Dissociation Energy De==DisE(I) & its 
c percent error (D1tDE) for each I.
c----------------------------------------------------------------
         DisE(i) = 0.0
      DO 50 j = 1,I
         if (key .eq. 0) then
           DisE(i)  = DisE(i) + CON(j)*y**j 
         else
           DisE(i)  = DisE(i) + (-1.0)**(j-1) * CON(j)*y**j 
         endif
 50   CONTINUE
C

C     IF ( De .gt. 0.0 ) THEN
      IF ( Deinp .gt. 0.0 ) THEN
         if ( nn.eq.1 ) then
            DisEi(i)  = DisE(i)
           DltDEi(i)  = 100.0 * abs(DisE(i)-De)/De
         else if( nn.eq.2 ) then
             DisE0(i) = DisE(i)
             DltE0(i) = 100.0 * abs(DisE0(i)-De)/De
           if (DisE(i) .ne. 0.0) DisEc(i) = DisE(i) + W0
             DltDEc(i) = 100.0 * abs(DisEc(i)-De)/De
         endif
      ELSE
C
           DisE(2) = CON(1)*CON(1)/( 4.0*CON(2) )
             if (DisE(2) .lt. 0.0) DisE(2) = - DisE(2) 
C
         if ( nn.eq.1 .and. i .gt. 2) then
            DisEi(i)  = DisE(i)
           DltDEi(i)  = 100.0 * abs(DisE(i)-DisE(i-1))/DisE(i-1)
         else if( nn.eq.2 .and. i .gt. 2) then
             DisE0(i) = DisE(i)
             DltE0(i) = 100.0 * abs(DisE0(i)-DisE0(i-1))/DisE0(i-1)
           if (DisE(i) .ne. 0.0) DisEc(i) = DisE(i) + W0
             DltDEc(i) = 100.0 * abs(DisEc(i)-DisEc(i-1))/DisEc(i-1)
         endif
           if (DisE(i-1) .eq. 0.0) then
	       DltDEi(i) = 0.0
	        DltE0(i) = 0.0
	       DltDEc(i) = 0.0
           endif
      ENDIF
c
  60   CONTINUE       
c
 100  CONTINUE       
c-------------------------------------------------------       
c  Print out the value of Calculated Disociation Energy  
c-------------------------------------------------------       
      IF (kr .gt. 0) GOTO 150
c
        write( 6,300)
        write(35,300)
C     DO i=2,7      
      DO i=2,mcst
	    if (i .eq. 2 .and. DisEc(i) .eq. 0.0) DltDEc(i+1) = 0.0
        write( 6,350) i,DisEc(i),DltDEc(i),Numaxc(i),der0(i),der2(i)
        write(35,350) i,DisEc(i),DltDEc(i),Numaxc(i),der0(i),der2(i)
      ENDDO       
C       Decal = DisEc(7)
C       ErrDe = 100.0* abs( DisEc(7) - DisEc(6) )/DisEc(7)
c-
        Decal = DisEc(mcst)
        ErrDe = 100.0* abs( DisEc(mcst) - DisEc(mcst-1) )/DisEc(mcst)
c-
 150    if (kr .eq. 0) then
          write( 6,310) De
          write(35,310) De
        else
	    write(38,310) De
        endif
C     DO i=2,7      
      DO i=2,mcst
	    if (i .eq. 2 .and. DisEi(i) .eq. 0.0) DltDEi(i+1) = 0.0
       if (kr .eq. 0) then
        write( 6,350) i,DisEi(i),DltDEi(i),Numaxi(i),der1(i),der3(i)
        write(35,350) i,DisEi(i),DltDEi(i),Numaxi(i),der1(i),der3(i)
       else
        write(38,350) i,DisEi(i),dltDEi(i),Numaxi(i),der1(i),der3(i)
       endif
      ENDDO       
C         Decal = DisEi(7)
C         ErrDe = 100.0* abs( DisEi(7) - DisEi(6) )/DisEi(7)
c-
          Decal = DisEi(mcst)
          ErrDe = 100.0* abs( DisEi(mcst) - DisEi(mcst-1) )/DisEi(mcst)
        IF (kr .gt. 0) GOTO 170
c-
        write( 6,320)
        write(35,320)
C     DO i=2,7      
      DO i=2,mcst
	    if (i .eq. 2 .and. DisE0(i) .eq. 0.0) DltE0(i+1) = 0.0
        write( 6,350) i, DisE0(i),DltE0(i),Numaxc(i),der0(i),der2(i)
        write(35,350) i, DisE0(i),DltE0(i),Numaxc(i),der0(i),der2(i)
      ENDDO
C         Decal = DisE0(7)
C         ErrDe = 100.0* abs( DisE0(7) - DisE0(6) )/DisE0(7)
c-
          Decal = DisE0(mcst)
          ErrDe = 100.0* abs( DisE0(mcst) - DisE0(mcst-1) )/DisE0(mcst)
 170    if (Numaxc(mcst) .gt. nv) nvp = Numaxc(mcst) + 1
c--------------------------------------------------------------------
 280  format(//5x,"NO REAL root for x=v+1/2 between  'star'  ",
     #"and  'endv', "/5x,"Please change the END value ! ",/
     #/15x,'v_start',8x,'v_end',/15x,f7.3,7x,f9.3,) 
 300  format(///12x,'**** Dissociation  Energy   WITH  W0 **** ',/
     #/2x,'n',5x,'De(cal)',8x,'ERROR_c%',4x,'v_max',5x,'dEv/dv',
     #8x,"dEv'/dv' ",/)
 310  format(//5x,'n is the number of vibrational constants used',
     #' to calc. De ',/
     #/8x,"Since De =/= De(w0,We+we0,WeXe,WeYe,...) is NOT a ",
     #/8x,"function of w0 & we0, so (We,WeXe,WeYe,WeZe, ...)",
     #/8x,"from (w0,We+we0,WeXe,WeYe,..) carrys ERROR in De_n !",/
     #/8x,'ERROR_n% = 100.0*abs( De_n - De )/De     FOR De=/=0 ',
     #/8x,'ERROR_n% = 100.0*|De_n - De_n-1|/De_n-1  FOR De===0 ',/
     #/11x,'For the INPUT  De =',1F15.9,3x,'a.u.'/
c    #/2x,'n',5x,'De(cal,n)',6x,'ERROR_n%',4x,'v_max',
     #/2x,'n',8x,'De_n ',8x,'ERROR_n%',4x,'v_max',
     #5x,'dEv/dv',8x,"dEv'/dv' ",/)
 320  format(///12x,'**** Dissociation Energy withOUT  W0 **** ',/
     #/2x,'n',5x,'De(cal)',8x,'ERROR_c%',4x,'v_max',5x,'dEv/dv',
     #8x,"dEv'/dv' ",/)
 350  format(1x,I2,2x,1pE15.9,2x,1pE11.5,2x,I4,3x,1pE11.4,4x,1pE11.4)    
c--------------------------------------------------------------------
         return
       END 
C===
      subroutine  calVIBconst(EV,mds)
      IMPLICIT  REAL*8(A-H,O-Z)
c--------------------------------------------------------------------
      logical   check
      DIMENSION ax(96,96),ex(30,96),fc(30,96),sv(30),EV(mds)
      DIMENSION xadev(30),eb(200),xn(40),ffn(40),wvib(40)
      DIMENSION xvar0(30),xsdev(30),xskew(30),xcurt(30),evi(200)
      DIMENSION eij(200,200),nj(200),evi0(200)
      common /LA1/ Deinp,Ner,mset,mcst,mtd,mdp,mv,nvs(20,80)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectrb/ Be,alphae,gamae,Der,betae,bryd,wvib0(40)
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /Eswitch/ aye,aze,ate,ase,are
      common /Eswitc3/ bxe,bye,bze,bte,bse,bre
      common /Eswitc4/ cxe,cye,cze,cte,cse,cre
      common /hnk/ hs2,Nryd,ksh
      common /fncom/ xmean(30),gf(30)
      common /spectr0/ w0,we0,wei,wex,wey,wez,wet,wes,wer,v0,r0
      common /spectr1/ Bee,ale,gae,eta3,eta4,eta5,eta6,eta7
      common /spectr2/ Dee,bete,xsi2,xsi3,xsi4,xsi5,xsi6,xsi7
      common /spectr3/ wen,wxe,wye,wze,wte,wse,wre
      common /spectr4/ w08,w09,w10,w11,w12,w13,w14,w15,w16
C-------------------------------------------------------------    
C   Generate the coefficients matrix ax(nv,nv) according to 
C Ev = We(v+1/2) - WeXe(v+1/2)**2 + WeYe(v+1/2)**3 + ...
C for (nv = nve+1) states :
c         v = 0, 1, 2, 3, 4, 5, ..., 10, ..., mv
C  mv =<  96 ;  ==> v_max = mv - 1 = 95.
C-------------------------------------------------------------  
          aucm = 219474.6306d0
        aJtoau = 0.229371d0
        aotoA0 = 0.5291772490d0
C-
        do n = 1, 200
            eb(n) = 0.0
          if (n .le. nvd) eb(n) = EV(n)
	  enddo
C-
      do 60 i=1, mv
          iv = i-1
        do k=1, mv
          if (nw0 .ge. 1) then
            ax(i,k) = (iv*1.0 + 0.5d0)**(k-1)
          else
            ax(i,k) = (iv*1.0 + 0.5d0)**k
          endif
        enddo
 60   continue
C-
        if (Deinp .eq. 0.0 .and. De .lt. 0.0) De = - De
C-------------------------------------------------------------    
C   Solve linear equation    ax * X = eb   to find the vector
C X = (W0, We+We0, WeXe, WeYe, WeZe, WeTe, WeSe, WeRe, ...)
C   When returns,  ex contains mset eb (solution vector X).
C-------------------------------------------------------------    
        call linersolv(ax,96,ex,30,96,200,eb,1)
C-------------------------------------------------------------    
C Cal. statistical data for every element of vector X (nw0=2):
C  X = (W0, We+We0, WeXe, WeYe, WeZe, WeTe, WeSe, WeRe, W08, ... W15)
C  i =   1,   2,      3,    4,    5,    6,   7,    8,    9,  ...  15
C                 We from INPUT;  We0 from calculations.
C     OR (for  nw0 = 1) :
C  X = (W0, We', WeXe, WeYe, WeZe, WeTe, WeSe, WeRe, W08, ... W15)
C  i =   1,  2,    3,    4,    5,    6,    7,    8,   8,  ...  15
C                 We'= We+We0  from solving  A*X = E.
C     OR (for  nw0 = 0) :
C  X = (We, WeXe, WeYe, WeZe, WeTe, WeSe, WeRe, W08, ... W15)
C  i =   1,   2,    3,    4,    5,    6,    7,    8, ...  15
C ex(k,i) -- the ith element generated from the kth data set.
C
C   Given an array of sv(n), moment returns its mean ave,
C average deviation adev, standard deviation sdev, variance
C var, skewness skew, and kurtosis curt.
C-------------------------------------------------------------    
      do 70 i=1, mcst
            do k=1, mset
               sv(k) = ex(k,i)
            enddo
          call moment(sv,mset,30,ave,adev,sdev,var,skew,curt)
            xmean(i) = ave
            xadev(i) = adev
            xsdev(i) = sdev
            xvar0(i) = var
            xskew(i) = skew
            xcurt(i) = curt
 70   continue
C---
        if (nw0 .ge. 1) write(38,330)
        if (nw0 .lt. 1) write(38,332)
      do 80 i=1, mcst
        write(38,340) i, xmean(i), xadev(i), xsdev(i), xvar0(i),
     # xskew(i), xcurt(i)
 80   continue
C---
          if (nlin .eq. 1 .and. nw0 .ge. 1) write(31,422)
          if (nlin .eq. 1 .and. nw0 .lt. 1) write(31,622)
          if (nlin .eq. 2 .and. nw0 .ge. 1) write(31,424)
          if (nlin .eq. 2 .and. nw0 .lt. 1) write(31,624)
        write(38,350) 
      do 90 i=1, mcst
        write(38,360) i, xmean(i), xmean(i)*aucm
        write(31,360) i, xmean(i), xmean(i)*aucm 
          if (i .eq. 8) write(31,*)
          if (i .eq. 8) write(38,*)
 90   continue
           w0 = 0.0
          we0 = 0.0
        if (nw0 .ge. 1)  w0 = xmean(1)
C-------------------------------------------------------------    
C  Cal. Ev using the AVERAGE values of vib. constants in xmean
C-------------------------------------------------------------    
          call  vibev(nvd,0,kn,200,evi)
C-
        if (nw0 .ge. 1)  write(41,450) De, bryd 
        if (nw0 .lt. 1)  write(41,650) De, bryd
          write(38,370) 
	      esum = 0.0
	      difs = 0.0
        do i=1, nvd
	      vdi0 = abs( eb(i) - evi(i) )
	      vdif = 100.0 * vdi0/eb(i)
	      esum = esum + vdif
	      difs = difs + vdi0
	      edif = evi(i) - evi(i-1) 
          write(38,360) i-1, eb(i), evi(i), vdif, edif
        enddo
          write(38,392)     esum/nvd
C----------------------------------
          write(40,432) De, difs/nvd
          write(41,434)  0, difs/nvd, 0.0, 0.0, 0.0
C-
          write(38,394) 
		   mk=0
        do i=1, nvd+kn
          write(38,360) i-1, eb(i), evi(i)
            if (i .eq. nvd) write(38,*)
            if (i .eq. nvd+5 .and. evi(i) .gt. evi(nvd) ) mk=1
        enddo
		if (mk .gt. 0) write(38,396)
C-------------------------------------------------------------    
C   ??? CHANGE the signs of wye, wte & wre by COMPARING with
C the ones obtained from using PERTURBATION formula.
C   Since    wen = We + We0,   so  We0 = wen - We .
C-------------------------------------------------------------    
           if (nw0 .ge. 1) then
               if ( w0 .eq. 0.0)  w0 = xmean(1)
              w0a = xmean(1)
              wen = xmean(2)
              wxe = xmean(3)
               if (we0 .eq. 0.0) we0 =  wen - We 
C-
c             wye =-xmean(4)
c             wze = xmean(5)
c             wte =-xmean(6)
c             wse = xmean(7)
c             wre =-xmean(8)
c             w09 = xmean(9)
c             w10 =-xmean(10)
c             w11 = xmean(11)
c             w12 =-xmean(12)
c             w13 = xmean(13)
c             w14 =-xmean(14)
c             w15 = xmean(15)
c             w16 =-xmean(16)
C-
              wye = xmean(4)
              wze = xmean(5)
              wte = xmean(6)
              wse = xmean(7)
              wre = xmean(8)
C-
              w08 = 0.0
              w09 = xmean(9)
              w10 = xmean(10)
              w11 = xmean(11)
              w12 = xmean(12)
              w13 = xmean(13)
              w14 = xmean(14)
              w15 = xmean(15)
              w16 = xmean(16)
           else
              wen = xmean(1)
              wxe = xmean(2)
              wye = xmean(3)
              wze = xmean(4)
              wte = xmean(5)
              wse = xmean(6)
              wre = xmean(7)
              w08 = xmean(8)
              w09 = xmean(9)
              w10 = xmean(10)
              w11 = xmean(11)
              w12 = xmean(12)
              w13 = xmean(13)
              w14 = xmean(14)
              w15 = xmean(15)
	     endif
C-
          call  FindnumDe(1,1,nvd,nvp,Decal,ErrDe)
                write(38,407) De, Decal, ErrDe, De-Decal
C-------------------------------------------------------------    
C   For (INPUT value of) De=0.0,  calculate De from
C         De == De( wen, wxe; wye)   
C   and generate initial GUESSES for FORCE constants 
C   fn's  in array "gf" using "guessFn".
C       De is used in "fmorse", "fryd" of "guessFn".
C-------------------------------------------------------------    
      IF ( Deinp .eq. 0.0 ) THEN
             write( 6,352) 
             write(38,352) 
          De2a = wen*wen/(4.0*wxe)
          De2b = We * We/(4.0*wxe)
           cubcal1 = wxe**2 - 3.0*wen*wye
           cubcal2 = wxe**2 - 3.0*We *wye
              e3a1 = 2.0*( dsqrt(cubcal1) )**3
              e3a2 = 2.0*( dsqrt(cubcal2) )**3
              e3b1 = wxe*( 2.0*wxe**2 - 9.0*wen*wye)
              e3b2 = wxe*( 2.0*wxe**2 - 9.0*We *wye)
          De3a = ( e3a1 - e3b1 )/(27.0*wye**2) 
          De3b = ( e3a2 - e3b2 )/(27.0*wye**2)
c
            ide = 1
		   if (De2a .lt. 0.0) De2a = - De2a
		   if (De2b .lt. 0.0) De2b = - De2b
		   if (De3a .lt. 0.0) De3a = - De3a
		   if (De3b .lt. 0.0) De3b = - De3b
          write( 6,342) De2a,De2b,De3a,De3b
          write(38,342) De2a,De2b,De3a,De3b
	       GOTO 900
c-
      ENDIF
C-------------------------------------------------------------    
 95       if (De .lt. 0.0) De = - De
C-------------------------------------------------------------    
C   Since the WeXe == wxe obtained from solving A*X=E  is 
C < 0.0 & this makes De < 0.0, so one HAS to change sign ! 
C
C   Using "guessFn" to generate INITIAL guess 'gf' for 
C FORCE constants  fn's.   Here  wen = We + we0 .
C-------------------------------------------------------------    
              call  guessFn(wen)
C-------------------------------------------------------------    
        IF (nlin .gt. 0) THEN
c
              if (nlin .eq. 1) write(33,412) bryd,De,wen
              if (nlin .eq. 2) write(33,414) bryd,De,wen
c
               gf2 = amu*wen*wen
               fa2 = gf2/(aJtoau * aotoA0**j)
c
              write(38,386) bryd, wen
		  write(38,400) 2, gf2, fa2
		  write(33,400) 2, gf2, fa2
            do i=1, 6
C           do i=1, mcst-1
                  j = i + 2
              xn(i) = gf(j)
                    faj = gf(j)/(aJtoau * aotoA0**j)
		      write(38,400) j, gf(j), faj
		      write(33,400) j, gf(j), faj
            enddo 
C=============================================================
C   Cal. FORCE constants fn's from the vib. constants in 
c xmean using Sun-Hou-Feng formulae by solving NON-linear
c MULTI-dimensional (6-D) equations (in subroutine) FUNCV
c using code "broydn"  for a given initial guess X(1:n)==xn 
c which is from early calculations made by "guessFn".
c 
c   "broydn" => "fdjac" => "funcv" --> 'fvec'(== Fn's). 
c 
C         call  broydn(xn,40,mcst-1,check)
c 
c   Since in the Sun-Hou-Feng's PERTURBATION formula, 6 
C FORCE constants (f3, f4, f5, f6, f7, f8) are expressed as
C the functions of 6 (+ 1) VIBrational spectrum constants
C      (We;  WeXe, WeYe, WeZe, WeTe, WeSe, WeRe)
C f9, f10, ... etc. are NOT defined in the formula, therefore,
C if mcst-1 > 6, "broydn" & "fvec" will solve a set of (mcst-1)
C dimension equations, and you'll get WRONG results ! 
C---------------------------------------------------------------    
          call  broydn(xn,40,6,check)
c
            do i=1, 6
C           do i=1, mcst-1
 	        ffn(i+2) = xn(i)
            enddo 
                write(38,398) 
		    write(33,630)
              ffn(2) = amu*wen*wen 
            do i=2, 8
C           do i=2, mcst+1
                faj = ffn(i)/(aJtoau * aotoA0**i)
              write(38,400) i, ffn(i), faj 
              write(33,400) i, ffn(i), faj 
                if (i .eq. 8) write(33,*)
                if (i .eq. 8) write(38,*)
            enddo 
c---
		  gf(2) = ffn(2)
                write(33,640) 
                write(38,640) 
                  sumdif = 0.0
C           do i=2, mcst+1
            do i=2, 8
                fdif = gf(i) - ffn(i)
                  sumdif = sumdif + abs( fdif )
                faf =  gf(i)/(aJtoau * aotoA0**i)
                faj = ffn(i)/(aJtoau * aotoA0**i)
              write(38,400) i, faf, faj, fdif
              write(33,400) i, faf, faj, fdif
C             write(38,400) i, gf(i), ffn(i), fdif 
C             write(33,400) i, gf(i), ffn(i), fdif
            enddo 
              write(38,645) sumdif/7
              write(33,645) sumdif/7
C-------------------------------------------------------------    
C             write(33,645) sumdif/mcst
C   Calculate vibrational spectrum constants 'wvib' using 
C FORCE constants 'ffn' generated above & PERTURBATION Eqs.
C-------------------------------------------------------------    
          call  vibconst(wen,ffn,wvib,40)
c
                 w0 = wvib(30)
                we0 = wvib(31)
c
              write(38,410) 
            do i=1, mcst
                vdif = 100.0 * abs( xmean(i) - wvib(i) )/xmean(i)
              write(38,360) i, xmean(i), wvib(i), vdif
C             write(31,360) i, wvib(i)
                xmean(i) = wvib(i)
            enddo
              write(38,550) w0a, w0, we0  
          if (nw0 .gt. 0 .and. w0a .ne. 0.0) w0 = w0a
C-------------------------------------------------------------    
C   Calculate vibrational ENERGIES 'evi' using vibrational
C constants 'wvib=xmean' generated above.
C-------------------------------------------------------------    
c             cxe = 1.0
c             cye = aye
c             cze = aze
c             cte = ate
c             cse = ase
c             cre = are
c---
              cxe = bxe
              cye = bye
              cze = bze
              cte = bte
              cse = bse
              cre = bre
c
            call  vibev(nvd,0,kn,200,evi)
c
              write(38,420) 
	          esum = 0.0
	          difs = 0.0
            do i=1, nvd
	        vdi0 = abs( eb(i) - evi(i) )
	        vdif = 100.0 * vdi0/eb(i)
	          difs = difs + vdi0
	          esum = esum + vdif
	          edif = evi(i) - evi(i-1) 
              write(38,360) i-1, eb(i), evi(i), vdif, edif
            enddo
              write(38,430)     esum/nvd
C             write(40,432) De, difs/nvd
C-
              write(38,394) 
                mk = 0
            do i=1, nvd+kn
              write(38,360) i-1, eb(i), evi(i)
                if (mk .eq. 0 .and. evi(i) .gt. evi(i+1) ) then
                  if ( evi(i) .gt. 2.0*evi(nvd) ) then
                      mk = 2
                  else
                      mk = 1
                    maxv = i-1
                     Dem = evi(i)
                     write(38,*)
                  endif
                endif
C                   if (mk .eq. 0 .and. i .eq. nvd) write(38,*)
                  edif = evi(i) - evi(i-1)
                  edi1 = evi(i+1) - evi(i)
                if (mk .eq. 0 .and. i .ge. 2) then
 	            if ( edi1 .gt. 2.0*edif ) then
		          mk=2
                    write(38,*)
                  endif
                endif
            enddo
                if (mk .gt. 1) write(38,396)
            if (mk .le. 1) write(38,405) maxv,Dem,Deinp,Dem-Deinp 
C-------------------------------------------------------------    
C  VIB. consts "wvib" are obtained using PERTURBATION formula
C-------------------------------------------------------------    
           if (nw0 .ge. 1) then
               if ( w0 .eq. 0.0)  w0 = wvib(1)
              wen = wvib(2)
               if (we0 .eq. 0.0) we0 =  wen - We
              wxe = wvib(3)
              wye = wvib(4)
              wze = wvib(5)
              wte = wvib(6)
              wse = wvib(7)
              wre = wvib(8)
C-
              w08 = 0.0
              w09 = wvib(9)
              w10 = wvib(10)
              w11 = wvib(11)
              w12 = wvib(12)
              w13 = wvib(13)
              w14 = wvib(14)
              w15 = wvib(15)
              w16 = wvib(16)
           else
              wen = wvib(1)
              wxe = wvib(2)
              wye = wvib(3)
              wze = wvib(4)
              wte = wvib(5)
              wse = wvib(6)
              wre = wvib(7)
              w08 = wvib(8)
              w09 = wvib(9)
              w10 = wvib(10)
              w11 = wvib(11)
              w12 = wvib(12)
              w13 = wvib(13)
              w14 = wvib(14)
              w15 = wvib(15)
           endif
C-
            call  FindnumDe(1,0,nvd,nvp,Decal,ErrDe)
                write(38,407) De, Decal, ErrDe, De-Decal
C-
        ENDIF
C--------------------------------------------------------
c  Loop over mset data sets generated from INPUT Ev's
C========================================================
 100      ik1 = 50
      do 200 k=1, mset
            if (k .eq. 1) write(38,372)
          write(38,380) k

C Modified by Chunjun Shu 2004-05-16

          write(80,380) k
          write(33,382) k
          write(31,384) k
C--------------------------------------------------------
C  Print the k_th set vibrational constants
C--------------------------------------------------------
        do i=1, mcst
          write(38,360) i, ex(k,i), ex(k,i)*aucm

C Modified by Chunjun Shu 2004-05-16

          write(80,360) i, ex(k,i), ex(k,i)*aucm
          write(31,360) i, ex(k,i), ex(k,i)*aucm
            if (i .eq. 7) write(38,*)
            if (i .eq. 7) write(80,*)
            if (i .eq. 7) write(31,*)
              xmean(i) = ex(k,i)
        enddo
            if (nw0 .gt. 0)  then
		   w0 = xmean(1)
		  w0a = w0
              we0 = xmean(2) - We
            else
		   w0 = 0.0
              we0 = xmean(1) - We
            endif
C--------------------------------------------------------
C  Generate & Compare Vibrational energies using ex(k,i)
C--------------------------------------------------------
            cxe = bxe
            cye = bye
            cze = bze
            cte = bte
            cse = bse
            cre = bre
C-
          call  vibev(nvd,0,kn,200,evi)
C-----------------------------------------------
C  Print & compare energies :
C-----------------------------------------------
          call  printev(nvd,kn,200,mk,maxv,esum,Dem,eb,evi)

C--------------------------------
C           write(80, ) is modified by Chunjun Shu
C---------------------------------

              aucm = 219474.6306d0

              sDem = Dem * aucm 
              sDeinp = Deinp * aucm 
            if (mk .gt. 1) write(38,396)
            if (mk .gt. 1) write(80,396)

            if (mk .le. 1) write(38,405) maxv,Dem,Deinp,Dem-Deinp 
                           write(38,392)    esum/nvd
                             esum0 = esum

C           if (mk .le. 1) write(80,406) maxv,Dem,Deinp,
C    $                           100 * (Dem-Deinp)/Deinp 
            if (mk .le. 1) write(80,406) maxv,sDem,sDeinp,
     $                           100 * (Dem-Deinp)/Deinp 
                           write(80,392)    esum/nvd
                             esum0 = esum
C              write(41,434)    k,esum/nvd,Dem,Deinp,Dem-Deinp
C              write(ik1+k,434) k,esum/nvd,Dem,Deinp,Dem-Deinp
C---
            write(38,375)
C---------------------------------
C Modified by Chunjun Shu
               open(unit=99,file="Dediff.tmp")
            if (mk .le. 1) then
               write(99,*)(Dem-Deinp)/Deinp
            else
            write(99,*) dabs(Dem-Deinp)
            endif
C Modified by Chunjun Shu
          open(unit=88,file="Ave.error.tmp")
          write(88,*)esum/nvd
          open(unit=77,file="DeEr.all")
          write(77,888)(Dem-Deinp)/Deinp,esum/nvd
888       format(2E15.8)
C----------------------------------
		  difsum = 0.0
        do i=1, nvd+kn
	    if (i .le. nvd) then
              evdif = eb(i)-evi(i)
            write(38,360) i-1, eb(i), evi(i), evdif
		  difsum = difsum + abs(evdif)
	    endif
              evi0(i) = evi(i)
        enddo
            write(38,377) difsum/nvd
              difsum0 = difsum
                  w0a = w0
                 we0a = we0
               write(41,434)    k,difsum/nvd,Dem,Deinp,Dem-Deinp
               write(ik1+k,434) k,difsum/nvd,Dem,Deinp,Dem-Deinp
C--------------------------------------------------------
C   Generate NUMerical DISSOCIATION energies De's 
C using the k_th set VIBrational constants from A*X=E.
C   wen = We + we0 ;      we0 = wen - We
C--------------------------------------------------------
           if (nw0 .ge. 1) then
               if ( w0 .eq. 0.0)  w0 = xmean(1)
              wen = xmean(2)
              wxe = xmean(3)
               if (we0 .eq. 0.0) we0 =  wen - We 
C-
c             wye =-xmean(4)
c             wze = xmean(5)
c             wte =-xmean(6)
c             wse = xmean(7)
c             wre =-xmean(8)
c             w09 = xmean(9)
c             w10 =-xmean(10)
c             w11 = xmean(11)
c             w12 =-xmean(12)
c             w13 = xmean(13)
c             w14 =-xmean(14)
c             w15 = xmean(15)
c             w16 =-xmean(16)
C-
c               write(38,92) wen,wxe,wye,wze,wte,wse,wre
C-
              wye = xmean(4)
              wze = xmean(5)
              wte = xmean(6)
              wse = xmean(7)
              wre = xmean(8)
C-
              w08 = 0.0
              w09 = xmean(9)
              w10 = xmean(10)
              w11 = xmean(11)
              w12 = xmean(12)
              w13 = xmean(13)
              w14 = xmean(14)
              w15 = xmean(15)
              w16 = xmean(16)
           else
              wen = xmean(1)
              wxe = xmean(2)
              wye = xmean(3)
              wze = xmean(4)
              wte = xmean(5)
              wse = xmean(6)
              wre = xmean(7)
              w08 = xmean(8)
              w09 = xmean(9)
              w10 = xmean(10)
              w11 = xmean(11)
              w12 = xmean(12)
              w13 = xmean(13)
              w14 = xmean(14)
              w15 = xmean(15)
           endif
C-
          if (k .eq. mdp) then
            do i=1, mcst+5
              wvib0(i) = xmean(i)
            enddo
              wvib0(30) = w0
              wvib0(31) = we0
          endif 
C-
            call  FindnumDe(1,1,nvd,nvp,Decal,ErrDe)
                write(38,407) De, Decal, ErrDe, De-Decal
C--------------------------------------------------------
C  Generate INITIAL guess for 'gf' for FORCE constants
C--------------------------------------------------------
            call  guessFn(wen)
C---
             ffn(2) = amu*wen*wen 
            fc(k,2) = ffn(2)
              gf(2) = ffn(2)
C---
C       IF (nlin .gt. 0) THEN
        IF (nlin .gt. 0 .and. Deinp .gt. 0.0) THEN
          do i=1, mcst-1
C         do i=1, mcst-1
            xn(i) = gf(i+2)
          enddo 
C-
            write(38,386) bryd, wen 
          do i=2, mcst+1
              faj = gf(i)/(aJtoau * aotoA0**i)
            write(38,400) i, gf(i), faj 
            write(33,400) i, gf(i), faj
              if (i .eq. 8) write(38,*)
              if (i .eq. 8) write(33,*)
          enddo 
C--------------------------------------------------------
C  Solve for FORCE constants ffn using gf as first guess
C           call  broydn(xn,40,mcst-1,check)
C--------------------------------------------------------
            call  broydn(xn,40,6,check)
C---
C         do i=1, mcst-1
          do i=1, 6
 	       ffn(i+2) = xn(i)
            fc(k,i+2) = xn(i)
          enddo 
            write(38,390) 
            write(33,630) 
C         do i=2, mcst+1
          do i=2, 8
              faj = ffn(i)/(aJtoau * aotoA0**i)
            write(38,400) i, ffn(i), faj 
            write(33,400) i, ffn(i), faj
              if (i .eq. 8) write(38,*)
              if (i .eq. 8) write(33,*)
          enddo 
c---
                write(33,640) 
                write(38,640) 
                  sumdif = 0.0
C           do i=2, mcst+1
            do i=2, 8
                fdif = gf(i) - ffn(i)
                  sumdif = sumdif + abs( fdif )
                faf =  gf(i)/(aJtoau * aotoA0**i)
                faj = ffn(i)/(aJtoau * aotoA0**i)
              write(38,400) i, faf, faj, fdif
              write(33,400) i, faf, faj, fdif
C             write(38,400) i, gf(i), ffn(i), fdif 
C             write(33,400) i, gf(i), ffn(i), fdif
            enddo 
              write(38,645) sumdif/7
              write(33,645) sumdif/7
C-------------------------------------------------------------    
C             write(33,645) sumdif/mcst
C   Calculate ROTational spectrum constants using 
C FORCE constants 'ffn' generated above.
C-------------------------------------------------------------    
            call  vjcoef(wen)
            call  convj(ffn,40,wen)
C-------------------------------------------------------------    
C   Print out calculated ROtational spectrum constants
C-------------------------------------------------------------    
            call  calEvj(1,mv,nj,0,200,evi,eij)
C-------------------------------------------------------------    
C   Calculate vibrational spectrum constants 'wvib' using 
C FORCE constants 'ffn' generated above.
C   w0a  is from solving  A*X = E  for nw0 > 0 .
C-------------------------------------------------------------    
            call  vibconst(wen,ffn,wvib,40)
C-------------------------------------------------------------    
              write(38,410) 
            do i=1, mcst
                vdif = 100.0 * abs( xmean(i) - wvib(i) )/xmean(i)
              write(38,360) i, xmean(i), wvib(i), vdif
C             write(31,360) i, wvib(i)
            enddo
C---
                 w0 = wvib(30)
		    we0 = wvib(31)
             if (nw0 .eq. 0) Wep = xmean(1)
             if (nw0 .gt. 0) Wep = xmean(2)
C-------------------------------------------------------------    
C   Print calculated  w0 = wvib(30); we0 = wvib(31);  and
C         We_solv = xmean(2)_solv - wvib(31) = We_pre - we0 .    
C-------------------------------------------------------------    
              write(38,362) wvib(30), wvib(31), Wep-wvib(31)
C-------------------------------------------------------------    
              write(38,550) w0a, w0, we0  
                if (nw0 .gt. 0) w0 = w0a
C-------------------------------------------------------------    
C   Calculate VIBrational energies using { We_solv, ... }.   
C-------------------------------------------------------------    
            call  vibev(nvd,1,kn,200,evi)
C-----------------------------------------------
C  Print & compare energies :
C-----------------------------------------------
          call  printev(nvd,kn,200,mk,maxv,esum,Dem,eb,evi)
                  write(38,392)    esum/nvd
                    esum1 = esum
C-------------------------------------------------------------    
C       We' = We_solv + we0  from solving   A*X = E .
C  For nw0 = 0 :
C    evi0(v) = We'*(v+1/2) - WeXe*(v+1/2)**2  + ...
C     evi(v) = w0 +  We'*(v+1/2) - WeXe*(v+1/2)**2  + ...
C       evhb = We_solv*(v+1/2) - WeXe*(v+1/2)**2  + ...
C
C  For nw0 > 0 :
C    evi0(v) = w0 + (We_solv +we0)*(v+1/2) - WeXe*(...) + ...
C     evi(v) = w0 +  We_solv*(v+1/2) - WeXe*(v+1/2)**2  + ...
C       evhb = (We_solv+we0)*(v+1/2) - WeXe*(v+1/2)**2  + ...
C
C  Compare ERROR% of energies :
C-------------------------------------------------------------    
        if (nw0 .eq. 0) write(38,520)
        if (nw0 .gt. 0) write(38,522)
	      difsum  = 0.0
	      difsum1 = 0.0
	      difsum2 = 0.0
              esum2 = 0.0
        do i=1, nvd
          if (nw0 .gt. 0) then 
            if (w0 .gt. 0.0) then 
              evhb = evi0(i) - w0a
C             evhb =  evi(i) - w0a
            else
              evhb = evi0(i) - cxe*w0a
C             evhb =  evi(i) - cxe*w0a
            endif
          else
              evhb = evi0(i) - we0
          endif
             evdif0 = eb(i)-evi0(i)
             evdif1 = eb(i)-evi(i)
             evdif2 = eb(i)-evhb
	      difsum  = difsum  + abs(evdif0)
	      difsum1 = difsum1 + abs(evdif1)
	      difsum2 = difsum2 + abs(evdif2)
	        esum2 = esum2 + 100.0*abs(evdif2)/eb(i) 
            write(38,360) i-1, eb(i), evdif0, evdif1, evdif2
	  enddo
            write(38,530) difsum/nvd, difsum1/nvd, difsum2/nvd
            write(38,540) esum0/nvd, esum1/nvd, esum2/nvd
C-------------------------------------------------------------    
c             cxe = 1.0
c             cye = aye
c             cze = aze
c             cte = ate
c             cse = ase
c             cre = are
c---
              cxe = bxe
              cye = bye
              cze = bze
              cte = bte
              cse = bse
              cre = bre
C--------------------------------------------------------
C   Generate & Compare Vibrational energies produced
c using ABOVE constants 'wvib'.
C           call  vibev(nvd,0,0,200,evi)
C--------------------------------------------------------
            do i=1, mcst
                xmean(i) = wvib(i)
            enddo
                xmean(20) = wvib(30)
                xmean(21) = wvib(31)
c---
            call  vibev(nvd,0,kn,200,evi)
c---
              write(38,420) 
	          esum = 0.0
	            mk = 0
            do i=1, nvd+kn
	            edif = evi(i) - evi(i-1) 
                if (mk .eq. 0 .and. edif .lt. 0.0) then
		      mk = 1
                  write(38,*) 
	          endif
	        if (i .le. nvd)  then
	          vdif = 100.0 * abs( eb(i) - evi(i) )/eb(i)
		    esum = esum + vdif
                write(38,360) i-1, eb(i), evi(i), vdif, edif
	        else
                write(38,360) i-1, eb(i), evi(i), 0.0, edif
	        endif
            enddo
c---
              write(38,430) esum/nvd
        ENDIF
C-
 200  continue
C=============================================================
C   Using the FORCE constants in "fc" obtained ABOVE from 
C the INPUT Ev's as the INItial guesses to CAlculate the 
C theoretical FORCE constants, VIBrational spectrum constants,
C and the VIBrational energies.
C------------------------------------------------------------- 
            write(38,500) 
      do 250 k=1, mset
C----------------------------------
C   Using "fc" as first guesses
C----------------------------------
            write(33,510) k
            write(38,510) k
C         do i=1,mcst-1
          do i=1, 6
             xn(i) = fc(k,i+2)
               faj = xn(i)/(aJtoau * aotoA0**i)
            write(33,400) i+2, xn(i), faj
            write(38,400) i+2, xn(i), faj
          enddo
C----------------------------------
C   Solve for FORCE constants ffn
C       call  broydn(xn,40,mcst-1,check)
C----------------------------------
        call  broydn(xn,40,6,check)
C-
              ffn(2) = fc(k,2)
C         do i=1,mcst-1
          do i=1, 6
             ffn(i+2) = xn(i)
          enddo
            write(33,390)
            write(38,390)
C         do i=2, mcst+1
          do i=2, 8
              faj = ffn(i)/(aJtoau * aotoA0**i)
            write(38,400) i, ffn(i), faj
            write(33,400) i, ffn(i), faj
              if (i .eq. 8) write(38,*)
              if (i .eq. 8) write(33,*)
          enddo
c---
                write(33,640) 
                write(38,640) 
                  sumdif = 0.0
C           do i=2, mcst+1
            do i=2, 8
               guesi = fc(k,i)
                fdif = guesi - ffn(i)
                  sumdif = sumdif + abs( fdif )
                faf = guesi/(aJtoau * aotoA0**i)
                faj = ffn(i)/(aJtoau * aotoA0**i)
              write(38,400) i, faf, faj, fdif
              write(33,400) i, faf, faj, fdif
C             write(38,400) i, guesi, ffn(i), fdif 
C             write(33,400) i, guesi, ffn(i), fdif
            enddo 
              write(38,645) sumdif/7
              write(33,645) sumdif/7
C------------------------------------------
C             write(33,645) sumdif/mcst
C   Cal. VIB. constants 'wvib' from 'ffn'
c      wei = We + we0 .
C------------------------------------------
           if (nw0 .ge. 1) then
             wei = ex(k,2)
           else
             wei = ex(k,1)
	     endif
C-
         call  vibconst(wei,ffn,wvib,40)
C-
            write(38,410)
          do i=1, mcst
              xmean(i) = ex(k,i)
            vdif = 100.0 * abs( xmean(i) - wvib(i) )/xmean(i)
              write(38,360) i, xmean(i), wvib(i), vdif
            xmean(i) = wvib(i)
          enddo
	       w0 = wvib(30)
	      we0 = wvib(31)
C-----------------------------------------
C   Cal. VIB. energies Ev's from 'wvib' 
C-----------------------------------------
          call  vibev(nvd,0,kn,200,evi)
c---
            write(38,420)
              esum = 0.0
                mk = 0
          do i=1, nvd+kn
              edif = evi(i) - evi(i-1)
                if (mk .eq. 0 .and. edif .lt. 0.0) then
		      mk = 1
                  write(38,*) 
	          endif
            if (i .le. nvd) then
              vdif = 100.0 * abs( eb(i) - evi(i) )/eb(i)
		  esum = esum + vdif
              write(38,360) i-1, eb(i), evi(i), vdif, edif
            else
              write(38,360) i-1, eb(i), evi(i), 0.0, edif
		endif
          enddo
            write(38,430) esum/nvd
C-
 250  continue
C=============================================================
 330  format(///6x,'***  Statistical calculations on molecular',
     #' constants  ***  : ',//
     #15x,' Ave -- mean (average) value of data; ',/
     #15x,'Adev -- average deviation of data; ',/
     #15x,'Sdev -- standard deviation of data; ',/
     #15x,' Var -- variance of data; ',/
     #15x,'Skew -- skewness of data; ',/
     #15x,'Curt -- kurtosis of data; ',//
     #9x,'{    W0,We+We0,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe,W08,...,W15 }',/
     #9x,'{ i = 1,   2,    3,   4,   5,   6,   7,   8,  9, ..., 15 }',/
     #//2x,'i',6x,'Ave',9x,'Adev',7x,'Sdev',7x,'Var',7x,
     #'Skew',8x,'Curt'/)
 332  format(///6x,'***  Statistical calculations on molecular',
     #' constants  ***  : ',//
     #15x,' Ave -- mean (average) value of data; ',/
     #15x,'Adev -- average deviation of data; ',/
     #15x,'Sdev -- standard deviation of data; ',/
     #15x,' Var -- variance of data; ',/
     #15x,'Skew -- skewness of data; ',/
     #15x,'Curt -- kurtosis of data; ',//
     #9x,'{    We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe,W08,...,W15 }',/
     #9x,'{ i = 1,  2,   3,   4,   5,   6,   7,   8,..., 15 }',/
     #//2x,'i',6x,'Ave',9x,'Adev',7x,'Sdev',7x,'Var',7x,
     #'Skew',8x,'Curt'/)
 340  format(1x,i2,1PE14.6,5(x,1PE10.3) )
 342  format(//1x,"ide =   1",12x,"2",12x,"3",12x,"4",
     #//2x,4(1pf12.7,x),/)
 350  format(//1x,'===  Average value of vibrational constants',
     #' ===',/
     #/1x,'These constants are obtained by SOLVING   A*X = E. ',
     #/1x,'NO ANAlytical formula for the constants for i > 7. ',
     #//3x,'i',5x,'Const(a.u.)',6x,'Const(cm-1)'/)
 352  format(//3x,"* Using  FOUR  De models for TEST calculations *",
     #//8x," wen, wxe, wye are from X = VIB_I : ",/
     #/11x,"ide = 1, De=De(wen,wxe) ",
     #/11x,"    = 2, De=De(We ,wxe) ",
     #/11x,"    = 3, De=De(wen,wxe,wye) ",
     #/11x,"    = 4, De=De(We ,wxe,wye) ",)
 360  format(2x,i2,4(1PE17.8) )
 362  format(/26x,'Const(a.u.)',//14x,'w0_cal =',1pe16.8,
     #/13x,'we0_cal =',1pe16.8,//13x,'We_solv =',1pe16.8,
     #' = We_pre - we0 ',/)
 365  format(//3x,'Initial guess for FORCE constants for : ',
     #//9x,'De_guess =',f15.10,//
     #2x,'(Obtained using De=De(We,WeXe) OR De=De(We,WeXe,WeYe) )',/
     #/2x,' n          f_(n; Har/ao**n)  ',/)
 370  format(//17x,'Comparing vibrational energies : ',/
     #/14x,'Einp(v) = Input vibrational energies; ',
     #/14x,'Ecal(v) = Calc. energies from above CONSTS. ',
     #/15x,'Error% = 100.0*| Einp(v)-Ecal(v)|/Einp(v)  ',
     #/14x,'Edif(v) = Ecal(v) - Ecal(v-1) ',/
     #/3x,'v',4x,'Einp(v; a.u.)',4x,'Ecal(v; a.u.)',
     #7x,'Error%',8x,'Edif(v; a.u.)',/)
 372  format(//1x,'===  Following are obtained using READED  Ev ===',)
 375  format(/17x,'CHECKING energy differences : ',/
     #/3x,'v',4x,'Einp(v; a.u.)',4x,'Ecal(v; a.u.)',
     #3x,'Einp(v)-Ecal(v) ',/)
 377  format(/2x,'Average energy difference, |Ave_dif| =',1pe16.8,/)
 380  format(//1x,'***  Vibrational constants from data set',
     #i3,'  ***',/6x,38(1H=),/
     #/1x,'These constants are obtained by SOLVING   A*X = E. ',
     #/1x,'NO ANAlytical formula for the constants for i > 7. ',
     #//3x,'i',5x,'Const(a.u.)',6x,'Const(cm-1)'/)
 382  format(//5x,'* Vib. FORCE constants from data set',i3,
     #//2x,' n     f_(n; Har/ao**n)     f_(n; aJ/A**n) ',/)
 384  format(//'* Conv. Vib. SPECTRUM const. from data set',i3,
     #//3x,'i',5x,'Const(a.u.)',6x,'Const(cm-1)',/)
C    #//6x,'i',5x,'Const(a.u.)',/)
 386  format(//6x,'INItial GUESSES of FORCE constants : ',
     #/4x,"( Obtained using 'guessFn', bryd, & We' ) ",/
     #/1x,'bryd =',f10.6," ;   We' = We_sol + we0 =",1pe13.5,/
     #/2x,' n     f_(n; Har/ao**n)     f_(n; aJ/A**n) ',/)
 390  format(//3x,'Vibrational FORCE constants for this set : ',
     #/3x,'(Obtained from above  GUESSES  & "broydn") ',/
     #/3x,'There are NO ANAlytical formula for n > 8 ',/
     #/2x,' n     f_(n; Har/ao**n)     f_(n; aJ/A**n) ',/)
 392  format(/5x,'Average error for this set, Ave_error% =',
     #1pe14.6,'%',//)
 394  format(3x,'v',4x,'Einp(v; a.u.)',4x,'Ecal(v; a.u.)',/)

 396  format(/3x,'ABOVE energy spectrum is NOT converged ! ',
     #/3x,'& the spectrum constants  are NOT good .',) 
 398  format(//3x,'Vibrational FORCE constants for this set : ',
     #/3x,'(Obtained from above  GUESSES  & "broydn") ',/
     #/3x,'There are NO ANAlytical formula for n > 8 ',/
     #/2x,' n     f_(n; Har/ao**n)     f_(n; aJ/A**n) ',/)
 400  format(2x,i2,3x,1PE18.10,2x,1PE18.10,2x,1PE14.6)
 402  format(2x,i2,3x,1PE28.20)
 
 405  format(/7x,'Energy spectrum MAY be  converged ! ',
     #/7x,'& the spectrum constants are GOOD .',/
     #/3x,'Maxi. vib. quantum number,  v_max =',i4,
     #/3x,'The theoretical  Ev_max { =< De } = ',1PE17.9,'  a.u.',
     #/3x,'The INPUT DISSOCIATION energy  De = ',1PE17.9,'  a.u.',
     #/3x,'                      Ev_max - De = ',1PE17.9,'  a.u.',) 

C----------------------------------------
C Modified by Chunjun Shu
 406  format(/7x,'Energy spectrum MAY be  converged ! ',
     #/7x,'& the spectrum constants are GOOD .',/
     #/3x,'Maxi. vib. quantum number,  v_max =',i4,
C    #/3x,'The theoretical  Ev_max { =< De } = ',1PE17.9,'  a.u.',
     #/3x,'The theoretical  Ev_max { =< De } = ',  F17.9,'  cm-1',
C    #/3x,'The INPUT DISSOCIATION energy  De = ',1PE17.9,'  a.u.',
     #/3x,'The INPUT DISSOCIATION energy  De = ',  F17.9,'  cm-1',
C    #/3x,'                      Ev_max - De = ',1PE17.9,'  %',) 
     #/3x,'                 (Ev_max - De)/De = ',1PE17.9,'  %',) 
C-----------------------------------------

 407  format(//8x,'ErrDe% = 100.0*| De_di0 |/De(cal,mcst)',
     # /8x,'De_di0 = De(cal,mcst) - De(cal,mcst-1) ',
     # /8x,'De_dif = De(try) - De(cal,mcst) ',
     #//6x,'De(try)       De(cal,mcst)        ErrDe% ',
     #10x,'De_dif ',/1PE16.8,1x,1PE16.8,1x,1PE16.8,1x,1PE16.8,/)
 410  format(//9x,'=*-  CHECK vibrational constants :  -*= ',
     #//3x,'  Error% = 100.0*| Pre_data - New_data |/Pre_data ',
     #//6x,' (New_data are calculated from FORCE consts) ',
     #//9x,' Pre_data         New_data ',
     # /3x,'i',4x,'Const(a.u.)',6x,'Const(a.u.)',9x,'Error%',/)
 412  format(//5x," FORCE constants from calculated Ev's : ",
     #/2x,'bryd =',f10.6,' ;  De =',f10.7,' ;  Wei =',1pe12.5,
     #//2x,' n     f_(n; Har/ao**n)     f_(n; aJ/A**n) ',/)
 414  format(//5x," FORCE constants from    INPUT   Ev's : ",
     #/2x,'bryd =',f10.6,' ;  De =',f10.7,' ;  Wei =',1pe12.5,
     #//2x,' n     f_(n; Har/ao**n)     f_(n; aJ/A**n) ',/)
 420  format(//18x,'CHECK vibrational energies : ',/
     #/14x,'Einp(v) = Input vibrational energies; ',
     #/14x,"Ecal(v) = Energies from CONSTS ( <=f_n's ). ",
     #/15x,'Error% = 100.0*| Einp(v)-Ecal(v)|/Einp(v)  ',
     #/14x,'Edif(v) = Ecal(v) - Ecal(v-1) ',/
     #/3x,'v',4x,'Einp(v; a.u.)',4x,'Ecal(v; a.u.)',
     #7x,'Error%',8x,'Edif(v; a.u.)',/)
 422  format(//" Converged VIB. SPECTRUM const. from calc. Ev's : "/
     #/1x,'{     W0, We+We0, WeXe, WeYe, WeZe, WeTe, WeSe, WeRe }',
     #/1x,'{ i =  1,    2,     3,    4,    5,    6,    7,    8  }',
     #//3x,'i',5x,'Const(a.u.)',6x,'Const(cm-1)',/)
 622  format(//" Converged VIB. SPECTRUM const. from calc. Ev's : "/
     #/1x,'{     We, WeXe, WeYe, WeZe, WeTe, WeSe, WeRe }',
     #/1x,'{ i =  1,   2,    3,    4,    5,    6,    7  }',
     #//3x,'i',5x,'Const(a.u.)',6x,'Const(cm-1)',/)
 424  format(//" Converged VIB. SPECTRUM const. from INPUT Ev's : "/
     #/1x,'{     W0, We+We0, WeXe, WeYe, WeZe, WeTe, WeSe, WeRe }',
     #/1x,'{ i =  1,    2,     3,    4,    5,    6,    7,    8  }',
     #//3x,'i',5x,'Const(a.u.)',6x,'Const(cm-1)',/)
 624  format(//" Converged VIB. SPECTRUM const. from INPUT Ev's : "/
     #/1x,'{     We, WeXe, WeYe, WeZe, WeTe, WeSe, WeRe }',
     #/1x,'{ i =  1,   2,    3,    4,    5,    6,    7  }',
     #//3x,'i',5x,'Const(a.u.)',6x,'Const(cm-1)',/)
 430  format(/5x,'Average error for this set, Ave_error% =',
     #1pe14.6,'%',/)
 630  format(/3x,'FORCE constants from "broydn" & above GUESS',
     #/2x,' n     f_(n; Har/ao**n)     f_(n; aJ/A**n) ',/)
 432  format(f12.7,6x,1PE17.9)
 434  format(2x,i2,1x,4(1PE14.6))
 440  format(//12x,' The following are obtained by solving : ', 
     #//8x,'       A * X = E_exp    where   X = VIB_I  ',
     #//8x," VIB_I --> De(try) & Fn_ini --> Fn's --> VIN_II ",
     #//8x,'De(try) = De(We,WeXe)_I   OR  De(We,WeXe,WeYe)_I ',
     #//8x,'De(cal) = De(We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe)_II ',
     #//9x,'De(inp)             De(try)              De(cal) ',
     #/3x,1PE17.9,3x,1PE17.9,3x,1PE17.9,/
     # /9x,' The OUTCOME from using  De(try) may be WRONG ! ',/)
C640  format(//9x,'Comparisons of FORCE constants { in  Har/ao**n } ',
 640  format(//9x,'Comparisons of FORCE constants { in  aJ/A**n } ',
     #/14x,'Dif_f(n) = f_(n; guess) - f_(n; calc.) ',//2x,
     #' n       f_(n; guess)        f_(n; calc.)        Dif_f(n)',/)
 645  format(/3x,"AVERAGE difference of | Dif_f(n) | is, ",
     #"sumdif =",1pe14.6,/)
 450  format(//14x,'The following are corresponds to ',
     #/9x,'De =',1pe16.8,' ;   bryd =',1pe13.5,/
     #/14x,'Best De should corresponds to a best   ',
     #/12x,'set of VIBrational constants & a VIB.  ',
     #/12x,'spectrum which has a SMALL |Ave_dif| ',
     #//15x,'|Ave_dif| = (SUM Dif_v)/(v_max+1) ',
     #//17x,'Dif_v = | Einp(v) - Ecal(v) | ',
C    #//12x,'Ave_error% = (SUM Error%_v)/(v_max+1) ',
C    #//12x,'Error%=100.0*| Einp(v)-Ecal(v)|/Einp(v) ',
C    #//15x,'De = De(We,WeXe,WeYe,WeZe,...) ',
     #/15x,'Ecal(v) = Ev(w0,We+we0,WeXe,WeYe,WeZe,...)  ',
     #//12x,'The k_th(set) is from the k_th solution ', 
     # /12x,'of the X(= W0,We+We0,WeXe,WeYe...) in  A*X = E. ',
     #//12x,'0_th set corresponds to the AVERAGE of  ',
     # /12x,'the constants X(= W0,We+We0,WeXe,WeYe,WeZe,...) ',/
C    #/2x,'k_th   Ave_err%    Ev_max(a.u.)   De(a.u)_inp ',
     #/2x,'k_th   |Ave_dif|   Ev_max(a.u.)   De(a.u)_inp ',
     #2x,'Ev_max - De ',/1x,'(set)',)
 650  format(//14x,'The following are corresponds to ',
     #/9x,'De =',1pe16.8,' ;   bryd =',1pe13.5,/
     #/14x,'Best De should corresponds to a best   ',
     #/12x,'set of VIBrational constants & a VIB.  ',
     #/12x,'spectrum which has a SMALL |Ave_dif| ',
     #//15x,'|Ave_dif| = (SUM Dif_v)/(v_max+1) ',
     #//17x,'Dif_v = | Einp(v) - Ecal(v) | ',
C    #/12x,'spectrum which has a SMALL Ave_error% ',
C    #//12x,'Ave_error% = (SUM Error%_v)/(v_max+1) ',
C    #//12x,'Error%=100.0*| Einp(v)-Ecal(v)|/Einp(v) ',
C    #//15x,'De = De(We,WeXe,WeYe,WeZe,...) ',
C    #/15x,'Ev = Ev(We,WeXe,WeYe,WeZe,...)  ',
     #/15x,'Ecal(v) = Ev(w0,We+we0,WeXe,WeYe,WeZe,...)  ',
     #//12x,'The k_th(set) is from the k_th solution ', 
     # /12x,'of the X( = We, WeXe, WeYe... ) in  A*X = E. ',
     #//12x,'0_th set corresponds to the AVERAGE of  ',
     # /12x,'the constants X( = We, WeXe, WeYe, WeZe,... ) ',/
     #/2x,'k_th   |Ave_dif|   Ev_max(a.u.)   De(a.u)_inp ',
     #2x,'Ev_max - De ',/1x,'(set)',)
 460  format(/6x,'Since De is required for INItial guess of FORCE ',
     #/6x,"constants Fn's, so NO Fn values as De_inp = 0.0",/)  
 500  format(//2x,55(1H=),
     #/2x,'  Using the FORCE constants in "fc" obtained ABOVE from ',
     #/2x,"the  INPUT   Ev's  as the  INItial guesses to CAlculate ", 
     #/2x,"the theoretical  FORCE constants,  VIBrational spectrum ", 
     #/2x,"constants, and the  VIBrational  energies. ",
     #//15x,"***  For REFERENCE ONLY  *** ",
     #/2x,55(1H-),/)
 510  format(//5x,'FIRST guess for FORCE consts from data set',i3,
     #/2x,"(Obtained from INPUT Ev's & PERturbation formula) ",/
     #/2x,'There are NO ANAlytical formula for n > 8 ',/
     #/2x,' n     f_(n; Har/ao**n)     f_(n; aJ/A**n) ',/)
 520  format(//17x,'Comparing energy differences & errors : ',/
     #/13x,'Einp(v) = Input vibrational energies . ',/
     #/14x,"   We' = We_solv + we0  from solving   A*X = E , ",
     #/14x,"     w0 adds POSITIVE contributions to E(v),  ",
     #/14x,"       Both w0 & we0 are IMPORTANT to E(v) ,  ",
     #/14x,"         we0 are from FORCE constants .  ",/
     #/14x,"Ev0(v) = We'*(v+1/2) - WeXe*(v+1/2)**2 + ... ",
     #/14x,"Ev1(v) = w0 + We'*(v+1/2) - WeXe*(...)**2 + ... ",
     #/14x,"Ev2(v) = We_solv*(v+1/2) - WeXe*(v+1/2)**2  + ... ",/
     #/14x,"The difference between Edif0 & Edif2 is about we0.",/
     #/14x,'Error% = 100.0*| Einp(v) - Evi(v)|/Einp(v)  ',
     #/20x,'Edifi(v) = Einp(v) - Evi(v) ',/
     #/3x,'v',4x,'Einp(v; a.u.)',3x,'Edif0(v;a.u.)',
     #4x,'Edif1(v;a.u.)',4x,'Edif2(v;a.u.)',/)
 522  format(//17x,'Comparing energy differences & errors : ',/
     #/13x,'Einp(v) = Input vibrational energies . ',/
     #/14x,"   We' = We_solv + we0  from solving   A*X = E , ",
     #/14x,"     w0 adds POSITIVE contributions to E(v),  ",
     #/14x,"       Both w0 & we0 are IMPORTANT to E(v) ,  ",
     #/14x,"         w0 & we0 are from FORCE constants .  ",/
     #/14x,"Ev0(v) = w0 + We'*(v+1/2) - WeXe*(...)**2 + ... ",
     #/14x,"Ev1(v) = w0 + We_solv*(v+1/2) - WeXe*(.)**2 + ... ",
     #/14x,"Ev2(v) = We'*(v+1/2) - WeXe*(v+1/2)**2  + ... ",/
C    #/14x,"Ev2(v) = We_solv*(v+1/2) - WeXe*(v+1/2)**2  + ... ",/
     #/14x,"The difference between Edif0 & Edif2 is about w0.",/
     #/14x,'Error% = 100.0*| Einp(v) - Evi(v)|/Einp(v)  ',
     #/20x,'Edifi(v) = Einp(v) - Evi(v) ',/
     #/3x,'v',4x,'Einp(v; a.u.)',3x,'Edif0(v;a.u.)',
     #4x,'Edif1(v;a.u.)',4x,'Edif2(v;a.u.)',/)
 530  format(/26x,'Average  Values  of  Differences & Errors',/
     #/11x,'Edifi_ave =',3(1pe16.8,x),)
 540  format(/10x,'Errori_ave =',3(1pe16.8,x),/)
 550  format(//14x,'CHECKING VIBrational constants  w0 : '/
     #/14x,'we0_calc --- from FORCE constants  ; ',
     #/15x,'w0_calc --- from FORCE constants  ; ',
     #/15x,'w0_solv --- from solving  A*X = E ; ',
     #/15x,'w0_solv =/= 0.0 ONLY for nw0 > 0 . ',/
     #/10x,'   w0_solv         w0_calc         we0_calc ',/
     #/7x,3(1pe16.8,x),/)
C
C-------------------------------------------------------------    
 900    return       
      END
C===
      subroutine vibev(nvi,ks,km,kv,evi)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
C-------------------------------------------------------------    
      dimension  evi(kv)
      common /spectr0/ w0,we0,wei,wex,wey,wez,wet,wes,wer,v0,r0
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectr4/ w08,w09,w10,w11,w12,w13,w14,w15,w16
      common /Eswitc4/ cxe,cye,cze,cte,cse,cre
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /fncom/ xmean(30),gf(30)
C-------------------------------------------------------------    
      if (nw0 .ge. 1) then
        if ( w0 .eq. 0.0)  w0 = xmean(1)
         Wei = xmean(2)
        wexi = xmean(3)
        weyi = xmean(4)
        wezi = xmean(5)
        weti = xmean(6)
        wesi = xmean(7)
        weri = xmean(8)
      else
         Wei = xmean(1)
        wexi = xmean(2)
        weyi = xmean(3)
        wezi = xmean(4)
        weti = xmean(5)
        wesi = xmean(6)
        weri = xmean(7)
	endif
          if (we0 .eq. 0.0) we0 = Wei - We
        if (ks .gt. 0 .and. w0 .eq. 0.0) then
            w0 = xmean(20)
           we0 = xmean(21)
	  endif
c---------------------------------------------
c  For ks > 0, We_solv = We' - we0_calc.  
c---------------------------------------------
      if (ks .eq. 0) then
C Modified by Chunjun Shu 2004-05-16
          write(38,10) w0,we0,Wei
          write(80,10) w0,we0,Wei
      else
              Wei = Wei - we0
          write(38,20) w0,we0,Wei
          write(80,20) w0,we0,Wei
	endif
c
          write(38,40) wexi,cxe,weyi,cye,wezi,cze,
     #weti,cte,wesi,cse,weri,cre
          write(80,40) wexi,cxe,weyi,cye,wezi,cze,
     #weti,cte,wesi,cse,weri,cre
c
  10   format(//13x,'---  In subroutine  "vibev"   --- ',/
     #/5x,' E(v) are calculated using following constants :',/ 
     #/5x,"{ For ks = 0 :  w0 = from fn's; or from  A*X = E ",
     #/5x,'               we0 = from FORCE constants;  OR ',
     #/5x,'               we0 = We_cal - We_inp = Wen - We  } ',/
     #/5x,"[ After soving A*X = E, you don't have fn's yet, ",
     #/5x,"  so,  w0 =/= 0 for nw0 > 0, w0 = 0 for nw0 = 0; ",
     #/5x,"      we0 has to be given by  we0 = Wen - We  .  ] ",/
     #/7x,'    w0 =',1PE18.10,
     #/7x,'   we0 =',1PE18.10,
     #/7x,'We_cal =',1PE18.10,)
  20   format(/13x,'==*  In subroutine  "vibev"   *== ',/
     #/9x,' E(v) are calculated using following constants :',/ 
     #/9x,"{ For ks > 0 :  w0 = from fn's; or from  A*X = E ",
     #/9x,"               we0 = from FORCE constants fn's ; ",
     #/19x,"We_sol = We' - weo_cal = Wei - we0      } ",/
     #/7x,'    w0 =',1PE18.10,
     #/7x,'   we0 =',1PE18.10,
     #/7x,'We_sol =',1PE18.10,)
  40   format(/7x,'  WeXe =',1PE18.10,5x,'cxe =',1pe10.2,
     #/7x,'  WeYe =',1PE18.10,5x,'cye =',1pe10.2,
     #/7x,'  WeZe =',1PE18.10,5x,'cze =',1pe10.2,
     #/7x,'  WeTe =',1PE18.10,5x,'cte =',1pe10.2,
     #/7x,'  WeSe =',1PE18.10,5x,'cse =',1pe10.2,
     #/7x,'  WeRe =',1PE18.10,5x,'cre =',1pe10.2,/)
c
      do 80 i=1,nvi+km
        evi(i) = 0.0d0
           bv = 1.0d0*( i - 1 )
          bv0 = bv + 0.5d0
         bv02 =  bv0*bv0
         bv03 = bv02*bv0
         bv04 = bv03*bv0
         bv05 = bv04*bv0
         bv06 = bv05*bv0
         bv07 = bv06*bv0
C-
         bv08 = bv07*bv0
         bv09 = bv08*bv0
         bv10 = bv09*bv0
         bv11 = bv10*bv0
         bv12 = bv11*bv0
         bv13 = bv12*bv0
         bv14 = bv13*bv0
         bv15 = bv14*bv0
C-
        if (nw0 .eq. 2) then
          if (w0 .lt. 0.0) then
            evi(i) =  cxe*w0 + (Wei + we0)*bv0 
          else
            evi(i) =  w0 + (Wei + we0)*bv0 
          endif
        elseif (nw0 .eq. 1) then
          if (w0 .lt. 0.0) then
            evi(i) =  cxe*w0 + Wei*bv0 
          else
            evi(i) =  w0 + Wei*bv0 
          endif
        elseif (nw0 .eq. 0) then
            evi(i) =  Wei*bv0 
        endif
C-
        if (wexi .lt. 0.0) then
          evi(i) =  evi(i) - cxe*wexi*bv02 
        else
          evi(i) =  evi(i) - wexi*bv02 
        endif
C-
        if (meny .eq. 0) then
          evi(i) =  evi(i) + cye*weyi*bv03 
          evi(i) =  evi(i) + cze*wezi*bv04 
          evi(i) =  evi(i) + cte*weti*bv05 
          evi(i) =  evi(i) + cse*wesi*bv06 
c
c--- Vibrational energies :
c
          evi(i) =  evi(i) + cre*weri*bv07
c-
c         evi(i) =  evi(i) + w08*bv08 + w09*bv09 + w10*bv10
c         evi(i) =  evi(i) + w11*bv11 + w12*bv12 + w13*bv13
c         evi(i) =  evi(i) + w14*bv14 + w15*bv15 
c--
        else
          evi(i) =  evi(i) + cye*weyi*bv03 
          evi(i) =  evi(i) - cze*wezi*bv04 
          evi(i) =  evi(i) + cte*weti*bv05 
          evi(i) =  evi(i) - cse*wesi*bv06 
          evi(i) =  evi(i) + cre*weri*bv07
c-
c         evi(i) =  evi(i) - w08*bv08 + w09*bv09 - w10*bv10
c         evi(i) =  evi(i) + w11*bv11 - w12*bv12 + w13*bv13
c         evi(i) =  evi(i) - w14*bv14 + w15*bv15 
        endif
c--
  80  continue
c--------------------------------------
        return
      END
C---
      subroutine printev(nvd,kn,mm,mk,maxv,esum,Dem,eb,evi)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      dimension  eb(mm),evi(mm)

C----------------------------
C Modified by Chunjun Shu 2004-05-03
      dimension  seb(200),sevi(200)
C------------------------------

C---
          write(38,100) 

C----------------------------
C Modified by Chunjun Shu 2004-05-03
          write(80,101) 
C----------------------------

	      esum = 0.0
	        mk = 0
        do i=1, nvd+kn
	      vdif = 0.0
	      edif = evi(i) - evi(i-1) 
	      edi1 = evi(i+1) - evi(i) 

C---------------------
C Modified by Chunjun Shu 2004-05-03

              aucm = 219474.6306d0
              seb(i) = eb(i) * aucm 
              sevi(i) = evi(i) * aucm 
              sedif = edif * aucm
C------------------------

	    if (i .le. nvd) then
	      vdif = 100.0 * abs( eb(i) - evi(i) )/eb(i)
	      esum = esum + vdif

C---------------------
C Modified by Chunjun Shu 2004-05-03

              write(38,150) i-1, eb(i), evi(i), vdif, edif
              write(80,151) i-1, seb(i), sevi(i), vdif, sedif
	     else
              write(38,150) i-1, eb(i), evi(i), 0.0, edif
              write(80,151) i-1, seb(i), sevi(i), 0.0, sedif
	    endif
C-------------------------

C             write(38,150) i-1, eb(i), evi(i), vdif, edif
C	    else

C             write(38,150) i-1, eb(i), evi(i), 0.0, edif
C	    endif
C-
            if (mk .eq. 0 .and. evi(i) .gt. evi(i+1) ) then
              if ( evi(i) .gt. 2.0*evi(nvd) ) then
                mk = 2
              else
                  mk = 1
                maxv = i-1
                 Dem = evi(i)
                 write(38,*)

                 write(80,*)

              endif
            endif
C                   if (mk .eq. 0 .and. i .eq. nvd) write(38,*)
              if (mk .eq. 0 .and. i .ge. 2) then
 	          if ( edi1 .gt. 2.0*edif ) then
		        mk=2
                  write(38,*)
                endif
              endif
        enddo
C---

 100  format(//17x,'Comparing vibrational energies : ',/
     #/14x,'Einp(v) = Input vibrational energies; ',
     #/14x,'Ecal(v) = Calc. energies from above CONSTS. ',
     #/15x,'Error% = 100.0*| Einp(v)-Ecal(v)|/Einp(v)  ',
     #/14x,'Edif(v) = Ecal(v) - Ecal(v-1) ',/
     #/3x,'v',4x,'Einp(v; a.u.)',4x,'Ecal(v; a.u.)',
     #7x,'Error%',8x,'Edif(v; a.u.)',/)
 101  format(//17x,'Comparing vibrational energies : ',/
     #/14x,'Einp(v) = Input vibrational energies; ',
     #/14x,'Ecal(v) = Calc. energies from above CONSTS. ',
     #/15x,'Error% = 100.0*| Einp(v)-Ecal(v)|/Einp(v)  ',
     #/14x,'Edif(v) = Ecal(v) - Ecal(v-1) ',/
     #/3x,'v',4x,'Einp(v; cm-1)',4x,'Ecal(v; cm-1)',
     #7x,'Error%',8x,'Edif(v; cm-1)',/)
 150  format(2x,i2,4(1PE17.8) )

C----------------------
C Modified by Chunjun Shu 2004-05-07
 151  format(2x,i2,F17.2,F17.3,2(1PE17.8) )

c--------------------------------------
        return
      END
C---
      subroutine linersolv(ax,ka,ex,ke,kf,mdd,bb,nk)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      dimension ax(ka,ka),ex(ke,kf),bb(mdd)
      dimension aa(96,96),ab(96,96),ac(96,96),bx(96,1)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /LA1/ Deinp,Ner,mset,mcst,mtd,mdp,mv,nvs(20,80)
C-------------------------------------------------
C Loop over (solve) the mset sets linear equa.
C       mset =< 20
C We INPUT only mcst nvs(i,j) for each i_th row.
C-------------------------------------------------
      do 30 j = 1, mset
        if (mtd .eq. 0) then
          m1 = nvs(j,1)
          m2 = m1 + mcst - 1
        else
          m1 = 1
          m2 = mcst
        endif
C--------------------------------------------------------
C For mtd=0 & i=m2 :  m2-m1+1 = m1+mcst-1 - m1+1 = mcst
C   The physical size of arrays :
C        aa(mcst,mcst),   bx(mcst,1)
C Define aa & bx arrays in the jth set.
C
C   The shifted energy  E(v) = Ev + De
C--------------------------------------------------------
        do 20 i = m1, m2
              i1 = i - m1 + 1
            if (mtd .eq. 0) then
              i2 = i
            else
              i2 = nvs(j,i)
            endif
          if (nk .eq. 0) bx(i1,1) = bb(i2) + De
          if (nk .gt. 0) bx(i1,1) = bb(i2)
C--------------------------------------------------------
C   Next line generates wrong  W0,We+We0,WeXe,WeYe, ... 
C for nk = 0 :     bx(i1,1) = bb(i2)
C    BUT,
C for nk = 1 :     bx(i1,1) = bb(i2)    is CORRECT !
C--------------------------------------------------------
          do k = 1, mcst
            aa(i1,k)=ax(i2,k)
            ab(i1,k)=ax(i2,k)
          enddo
 20     continue
            kt = 38
C-------------------------------------------------
C   Solve   aa*X=bx    for the jth set
C-------------------------------------------------
        call gaussj(aa,mcst,96,bx,1,1)
C-------------------------------------------------
C   When gaussj returns, aa = 1/ab , bx == X .
C Save the jth solution vector bx onto ex .
C-------------------------------------------------
            do k = 1, mcst
              ex(j,k) = bx(k,1)
            enddo
C-------------------------------------------------
C Print ab and its inverse matrix aa == 1/ab
C-------------------------------------------------
          write(kt,301) j
        call mprint(ab,mcst,mcst,96,96,kt,1)
        call mprint(aa,mcst,mcst,96,96,kt,2)
C-------------------------------------------------
C Check & print  if  ab *aa = ab * 1/ab = 1
C-------------------------------------------------
        call multp(mcst,mcst,mcst,aa,ab,ac,96,96,96)
        call mprint(ac,mcst,mcst,96,96,kt,3)
 30   continue
C-------------------------------------------------
C Print the solution vector array  ex
C-------------------------------------------------
        call mprint(ex,mset,mcst,30,96,kt,4)
C-------------------------------------------------
 301  format(//13x,'*** Solved A*X=b for data set i =',
     #i3,'  ***')
C-------------------------------------------------
c       return
      END
C===
      SUBROUTINE gaussj(a,n,np,b,m,mp)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
C============================================================================
C  GAUSSJ solves for the linear equations 
C                          A*X=b.  
C using the Gauss-Jordan elimination.    
C   Input:  a(n,n) is the coefficient matrix A.  b(n,m) is an inputting 
C           matrix containning the m right-hand side vectors.
C   Output: a(n,n) is the inverse matrix A-1 of A.  b(n,m) is the solution
c           vectors X; and b == a = A-1 if b is inputted as an unit matrix.
C============================================================================
        PARAMETER (NMAX=100)
        dimension  a(np,np),b(np,mp)
        INTEGER indxc(NMAX),indxr(NMAX),ipiv(NMAX)
C----------------------------------------------------------------------------
C
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
C 
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
C--------------------------------------------
        return
      END
C
      SUBROUTINE MULTP(N,L,M,A,B,C,N1,L1,M1)
	IMPLICIT REAL*8 (A-H,O-Z)
	dimension A(N1,L1),B(L1,M1),C(N1,M1)
C MATRIX MULTIPLICATION     C(N-M) = A(N-L) * B(L-M)
C
       DO 20 I=1,N
       DO 20 J=1,M
	C(I,J)=0.0d0 
       DO 20 K=1,L
	C(I,J) = C(I,J) + A(I,K)*B(K,J)
 20    CONTINUE
c------------------------------------------
	RETURN
	END
C
      subroutine mprint(g,m,n,ni,nj,kt,kg)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      dimension g(ni,nj)
c-------------------------------------------------------------------------
c  This subroutine is written to print a matrix g.  Weiguo Sun  06/25/1993
c On entry :
c   g(m,n)  -- The matrix to be printed.
c     m,n   -- Integers to specify the actual size of g.
c    ni,nj  -- Integers to specify the dimension of g.
c       kt  -- Index for printing unit.
c       kg  -- Switch to print the tittle (you may change it) of g. 
c              kg=0, no tittle is printed.
c-------------------------------------------------------------------------
      iout=kt
c     iout=6
c
        if (kg .eq. 0) go to 10
      if (kg .eq. 1) then
        write(iout,100)
      elseif (kg .eq. 2) then
        write(iout,110)
      elseif (kg .eq. 3) then
        write(iout,120)
      elseif (kg .eq. 4) then
        write(iout,130)
      endif
c
  10  do 70 k=1,n
        n1=k*5
      if(n1.ge.n) goto 75
        n2=n1-4
      write(iout,170) (j,j=n2,n1)
        write(iout,180)
C
      do 50 i=1,m
 50   write(iout,190) i,(g(i,j),j=n2,n1)
        write(iout,200)
 70   continue
C
 75     n2=n1-4
      write(iout,170) (j,j=n2,n)
      write(iout,180)
      do 80 i=1,m
 80   write(iout,190) i,(g(i,j),j=n2,n)
      write(iout,220)
c-------------------------------------------------------
 100  format(//23x,'The input matrix A :'/)
 110  format(//23x,'The inverse of A (=> 1/A) :'/)
 120  format(//23x,'The multiplication of A * (1/A) :'/)
 130  format(///23x,'The solution vector array  X : '/)
 170  format(4x,1hj,8x,i2,7(12x,i2))
 180  format(3x,1hi)
 190  format(2x,i2,2x,6(1pe14.6))
 200  format(1h*,33x,9hContinue./)
 220  format(1h!,35x,4hEND.)
c-------------------------------------------------------
        return
      end
C
      SUBROUTINE moment(data,n,nd,ave,adev,sdev,var,skew,curt)
C---------------------------------------------------------------------
C   Given an array of data(1:n), this routine returns its mean ave,
C average deviation adev, standard deviation sdev, variance var,
C skewness skew, and kurtosis curt.
C---------------------------------------------------------------------
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      dimension  data(nd)
C     REAL adev,ave,curt,sdev,skew,var,data(n)
C---------------------------------------------------------------------
      if(n.le.1)pause 'n must be at least 2 in moment'
      s=0.
      do 11 j=1,n
        s=s+data(j)
11    continue
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        adev=adev+abs(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
12    continue
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=sqrt(var)
      if(var.ne.0.)then
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3.
      else
        pause 'no skew or kurtosis when zero variance in moment'
      endif
        return
      END
C
C===
      SUBROUTINE broydn(x,n2,n,check)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
C----------------------------------------------------------------
c   To find the ROOTS of a user supplied MULTI-dimensional
c (n-D) equations (in subroutine) FUNCV for a given initial
c guess X(1:n).   "FUNCV" is called by "fdjac".
C   Given an initial guess X(1:n) for a ROOT in n dimensions,
C find the root by Broyden's method embedded in a globally 
C convergent strategy.  The vector of functions to be zeroed,
C called FVEC(1:n) in the routine below, is returned by a user
C supplied subroutine that MUST be called FUNCV and have the
C declaration subroutine  FUNCV(n,x,fvec). The subroutine FDJAC
C and the function FMIN are used.  
C   The output quantity CHECK is false on a normal return and 
C true if the routine has converged to a local minimum of the 
C function FMIN or if Broyden's can make no further progress.
C In this case try RESTARTING from a different initial guess.
C
C       PARAMETERS :
C   NP   -- is the maximum expected value of n;
C MAXITS -- is the maximum number of iterations;
C   EPS  -- is close to the machine precision;
C  TOLF  -- sets the convergence criterion on function values;
C TOLMIN -- sets the criterion for deciding whether spurious
C             convergence to a minimum of FMIN has occurred;
C  TOLX  -- is the convergence criterion on delta_x;
C STPMX  -- is the scaled maximum step length allowed in line 
C             searches.
C
C---   TOLF & EPS determines the quality and changes the 
C--- value of the FORCE constants  f_n's.
C
Cc    PARAMETER (n1=6, NP=40, MAXITS=200, EPS=1.e-7)
cc                                            *****
Cc    PARAMETER (n1=6, NP=40, MAXITS=200, EPS=1.e-9)
Cc    PARAMETER (TOLF=1.e-6, TOLMIN=1.e-8, TOLX=EPS, STPMX=100.0)
cc    PARAMETER (TOLF=1.e-4, TOLMIN=1.e-6, TOLX=EPS, STPMX=100.0)
cc                                  *****
C     dimension  x(n),c(NP),d(NP),fvcold(NP),g(NP),p(NP)
cc  
cc       Note :    n2 == NP  !
C----------------------------------------------------------------
      LOGICAL check,restrt,sing,skip
      PARAMETER (NP=40, MAXITS=200, EPS=1.e-7)
      PARAMETER (TOLF=1.e-4, TOLMIN=1.e-6, TOLX=EPS, STPMX=100.0)
      EXTERNAL fmin
      COMMON /newtv/ fvec(NP),nn
      dimension  x(n2),c(NP),d(NP),fvcold(NP),g(NP),p(NP)
      dimension  qt(NP,NP),r(NP,NP),s(NP),t(NP),w(NP),xold(NP)
C----------------------------------------------------------------
CU  USES fdjac,fmin,lnsrch,qrdcmp,qrupdt,rsolv
C----------------------------------------------------------------
        nn=n
      f=fmin(n2,n,x)
        test=0.
      do 11 i=1,n
        if(abs(fvec(i)).gt.test)test=abs(fvec(i))
11    continue
      if(test.lt..01*TOLF) return
      sum=0.
      do 12 i=1,n
        sum=sum+x(i)**2
12    continue
      stpmax=STPMX*max(sqrt(sum),float(n))
      restrt=.true.
      do 44 its=1,MAXITS
        if(restrt)then
          call fdjac(NP,n,x,fvec,r)
          call qrdcmp(r,n,NP,c,d,sing)
          if(sing) pause 'singular Jacobian in broydn'
          do 14 i=1,n
            do 13 j=1,n
              qt(i,j)=0.
13          continue
            qt(i,i)=1.
14        continue
          do 18 k=1,n-1
            if(c(k).ne.0.)then
              do 17 j=1,n
                sum=0.
                do 15 i=k,n
                  sum=sum+r(i,k)*qt(i,j)
15              continue
                sum=sum/c(k)
                do 16 i=k,n
                  qt(i,j)=qt(i,j)-sum*r(i,k)
16              continue
17            continue
            endif
18        continue
          do 21 i=1,n
            r(i,i)=d(i)
            do 19 j=1,i-1
              r(i,j)=0.
19          continue
21        continue
        else
          do 22 i=1,n
            s(i)=x(i)-xold(i)
22        continue
          do 24 i=1,n
            sum=0.
            do 23 j=i,n
              sum=sum+r(i,j)*s(j)
23          continue
            t(i)=sum
24        continue
          skip=.true.
          do 26 i=1,n
            sum=0.
            do 25 j=1,n
              sum=sum+qt(j,i)*t(j)
25          continue
            w(i)=fvec(i)-fvcold(i)-sum
            if(abs(w(i)).ge.EPS*(abs(fvec(i))+abs(fvcold(i))))then
              skip=.false.
            else
              w(i)=0.
            endif
26        continue
          if(.not.skip)then
            do 28 i=1,n
              sum=0.
              do 27 j=1,n
                sum=sum+qt(i,j)*w(j)
27            continue
              t(i)=sum
28          continue
            den=0.
            do 29 i=1,n
              den=den+s(i)**2
29          continue
            do 31 i=1,n
              s(i)=s(i)/den
31          continue
            call qrupdt(r,qt,n,NP,t,s)
            do 32 i=1,n
              if(r(i,i).eq.0.) pause 'r singular in broydn'
              d(i)=r(i,i)
32          continue
          endif
        endif
        do 34 i=1,n
          sum=0.
          do 33 j=1,n
            sum=sum+qt(i,j)*fvec(j)
33        continue
          g(i)=sum
34      continue
        do 36 i=n,1,-1
          sum=0.
          do 35 j=1,i
            sum=sum+r(j,i)*g(j)
35        continue
          g(i)=sum
36      continue
        do 37 i=1,n
          xold(i)=x(i)
          fvcold(i)=fvec(i)
37      continue
        fold=f
        do 39 i=1,n
          sum=0.
          do 38 j=1,n
            sum=sum+qt(i,j)*fvec(j)
38        continue
          p(i)=-sum
39      continue
        call rsolv(r,n,NP,d,p)
        call lnsrch(NP,n,xold,fold,g,p,x,f,stpmax,check)
C       call lnsrch(NP,n,xold,fold,g,p,x,f,stpmax,check,fmin)
        test=0.
        do 41 i=1,n
          if(abs(fvec(i)).gt.test)test=abs(fvec(i))
41      continue
        if(test.lt.TOLF)then
          check=.false.
          return
        endif
        if(check)then
          if(restrt)then
            return
          else
            test=0.
            den=max(f,.5*n)
            do 42 i=1,n
              temp=abs(g(i))*max(abs(x(i)),1.)/den
              if(temp.gt.test)test=temp
42          continue
            if(test.lt.TOLMIN)then
              return
            else
              restrt=.true.
            endif
          endif
        else
          restrt=.false.
          test=0.
          do 43 i=1,n
            temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
            if(temp.gt.test)test=temp
43        continue
          if(test.lt.TOLX)return
        endif
44    continue
        pause ' *  MAXITS exceeded in broydn  * '
      END
C
C================================================================
C
      SUBROUTINE fdjac(np,n,x,fvec,df)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
C---
C---   EPS determines the quality and changes the value of 
C--- the FORCE constants  f_n's.
C---
C     PARAMETER (NMAX=40,EPS=1.e-2)
      PARAMETER (NMAX=40,EPS=1.e-4)
C     PARAMETER (NMAX=40,EPS=1.e-8 )
      dimension  df(np,np),fvec(np),x(np),f(NMAX)
C----------------------------------------------------------------
      do 12 j=1,n
        temp=x(j)
        h=EPS*abs(temp)
        if(h.eq.0.)h=EPS
        x(j)=temp+h
        h=x(j)-temp
        call funcv(np,n,x,f)
C       call funcv(n,n,x,f)
        x(j)=temp
        do 11 i=1,n
          df(i,j)=(f(i)-fvec(i))/h
11      continue
12    continue
        return
      END
C
C================================================================
C
      FUNCTION fmin(n1,n,x)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      PARAMETER (NP=40)
C     PARAMETER (n2=6, NP=40)
      COMMON /newtv/ fvec(NP),nn
      dimension  x(n1)
c-------------------------------------------
c  Note :  n1 == NP
c-------------------------------------------
      call funcv(NP,n,x,fvec)
C     call funcv(n,n,x,fvec)
        sum=0.
      do 11 i=1,n
        sum=sum+fvec(i)**2
11    continue
      fmin=0.5*sum
        return
      END
C
C================================================================
c     SUBROUTINE lnsrch(np,n,xold,fold,g,p,x,f,stpmax,check,func1)
C
      SUBROUTINE lnsrch(np,n,xold,fold,g,p,x,f,stpmax,check)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      LOGICAL check
      PARAMETER (ALF=1.e-4, TOLX=1.e-7)
C     PARAMETER (ALF=1.e-9, TOLX=1.e-12)
      EXTERNAL fmin
      dimension  g(np),p(np),x(np),xold(np)
c-------------------------------------------------------------------
      check=.false.
        sum=0.
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=sqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      test=0.
      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        f=fmin(np,n,x)
C       f=fmin(n,x)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.1.)then
            tmplam=-slope/(2.*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.)then
              tmplam=-slope/(2.*b)
            else
              disc=b*b-3.*a*slope
              tmplam=(-b+sqrt(disc))/(3.*a)
            endif
            if(tmplam.gt..5*alam)tmplam=.5*alam
          endif
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=max(tmplam,.1*alam)
      goto 1
      END
C
C================================================================
C
      SUBROUTINE qrdcmp(a,n,np,c,d,sing)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      LOGICAL sing
      dimension  a(np,np),c(np),d(np)
c-------------------------------------------------------------------
      sing=.false.
      scale=0.
      do 17 k=1,n-1
        do 11 i=k,n
          scale=max(scale,abs(a(i,k)))
11      continue
        if(scale.eq.0.)then
          sing=.true.
          c(k)=0.
          d(k)=0.
        else
          do 12 i=k,n
            a(i,k)=a(i,k)/scale
12        continue
          sum=0.
          do 13 i=k,n
            sum=sum+a(i,k)**2
13        continue
          sigma=sign(sqrt(sum),a(k,k))
          a(k,k)=a(k,k)+sigma
          c(k)=sigma*a(k,k)
          d(k)=-scale*sigma
          do 16 j=k+1,n
            sum=0.
            do 14 i=k,n
              sum=sum+a(i,k)*a(i,j)
14          continue
            tau=sum/c(k)
            do 15 i=k,n
              a(i,j)=a(i,j)-tau*a(i,k)
15          continue
16        continue
        endif
17    continue
      d(n)=a(n,n)
      if(d(n).eq.0.)sing=.true.
      return
      END
C
C================================================================
C
      SUBROUTINE qrupdt(r,qt,n,np,u,v)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      dimension  r(np,np),qt(np,np),u(np),v(np)
c-------------------------------------------------------------------
      do 11 k=n,1,-1
        if(u(k).ne.0.)goto 1
11    continue
      k=1
1     do 12 i=k-1,1,-1
        call rotate(r,qt,n,np,i,u(i),-u(i+1))
        if(u(i).eq.0.)then
          u(i)=abs(u(i+1))
        else if(abs(u(i)).gt.abs(u(i+1)))then
          u(i)=abs(u(i))*sqrt(1.+(u(i+1)/u(i))**2)
        else
          u(i)=abs(u(i+1))*sqrt(1.+(u(i)/u(i+1))**2)
        endif
12    continue
      do 13 j=1,n
        r(1,j)=r(1,j)+u(1)*v(j)
13    continue
      do 14 i=1,k-1
        call rotate(r,qt,n,np,i,r(i,i),-r(i+1,i))
14    continue
      return
      END
C
C================================================================
C
      SUBROUTINE rsolv(a,n,np,d,b)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      dimension  a(np,np),b(np),d(np)
c-------------------------------------------------------------------
      b(n)=b(n)/d(n)
      do 12 i=n-1,1,-1
         sum=0.
         do 11 j=i+1,n
           sum=sum+a(i,j)*b(j)
 11      continue
         b(i)=(b(i)-sum)/d(i)
 12    continue
               return
      END
C===
      SUBROUTINE rotate(r,qt,n,np,i,a,b)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      dimension  r(np,np),qt(np,np)
c-------------------------------------------------------------------
      if(a.eq.0.)then
        c=0.
        s=sign(1.,b)
      else if(abs(a).gt.abs(b))then
        fact=b/a
        c=sign(1./sqrt(1.+fact**2),a)
        s=fact*c
      else
        fact=a/b
        s=sign(1./sqrt(1.+fact**2),b)
        c=fact*s
      endif
      do 11 j=i,n
        y=r(i,j)
        w=r(i+1,j)
        r(i,j)=c*y-s*w
        r(i+1,j)=s*y+c*w
11    continue
      do 12 j=1,n
        y=qt(i,j)
        w=qt(i+1,j)
        qt(i,j)=c*y-s*w
        qt(i+1,j)=s*y+c*w
12    continue
        return
      END
C
C===
      subroutine  getEvDe(nv,kw,lda,Lv,Ev)
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /spectr0/ w0,we0,wei,wex,wey,wez,wet,wes,wer,v0,r0
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /LA1/ Deinp,Ner,mset,mcst,mtd,mdp,mv,nvs(20,80)
      dimension  e1(200),Ew0(200),Ew(200),Ev(Lv)
c----------------------------------------------------------
c  nv = the number of vibrational states used.
c         (nv-1) - The HIGHEST vibrational state.
c==========================================================
       rydev = 13.60569809d0
        auev = 27.21139618d0
        aucm = 219474.6306d0
      rinert = amu*Re*Re
      Brigid = 1.0/(2.0*rinert)
c================================================================
c     do i=1,nv
      do i=1,nvd
c-- Read in the INPUT vibrational energies OR their differencies
        read(23,*) Ew0(i)
      enddo
c-
c     if (Ew0(nvd) .lt. Ew0(1)) then
      if (Ew0(nvd) .lt. Ew0(nvd-1)) then
c-- As Ew0(nvd) < Ew0(1), you readed Ev differencies;
c--   you need to find Ev themselves.
          Ew(1) = Ew0(1)
          Ew(2) = Ew(1) + Ew0(2)
        do i=3,nv 
          Ew(i) = Ew(i-1) + Ew0(i)
          if (i .gt. nvd)  Ew(i) = 0.0
          if (i .gt. nvd) Ew0(i) = 0.0
        enddo
      endif
c
C     if (Ew0(nvd) .gt. Ew0(1)) then
      if (Ew0(nvd) .gt. Ew0(nvd-1)) then
c-- As Ew0(nvd) > Ew0(1), you readed Ev ;
c--   you need to find Ev differencies .
          e1(1) = Ew0(1)
c       do i=2,nv
        do i=1,nvd
          e1(i) = Ew0(i) - Ew0(i-1)
        enddo
c  
        do i=1,nv
           Ew(i) = Ew0(i)
          Ew0(i) = e1(i)
           e1(i) = 0.0
          if (i .gt. nvd)  Ew(i) = 0.0
          if (i .gt. nvd) Ew0(i) = 0.0
        enddo
      endif
c
         if ( Deinp .eq. 0.0 ) Ev(1) = We*(0.0 + 0.5)
c
      if (Ew0(1) .gt. Ev(1)*1000.0d0) then
c-- Change the unit of Ew0 from cm-1 to a.u.
c       do i=1,nv
        do i=1,nvd
          Ew0(i) = Ew0(i)/aucm
           Ew(i) = Ew(i)/aucm
        enddo
      endif
c---------------------------------------------------------------- 
c  Print out the readed energies and their differencies
c---------------------------------------------------------------- 
        write( 6,500)
        write(38,500)
      do i=1,nvd
          edif = Ew0(i) - Ew0(i-1)
        write( 6,570) i-1, Ew(i), Ew0(i), edif
        write(38,570) i-1, Ew(i), Ew0(i), edif
          if ( Deinp .eq. 0.0 ) Ev(i) = Ew(i)
      enddo
        write( 6,510)
        write(38,510)
      do i=1,nvd
          edif = Ew0(i) - Ew0(i-1)
        write( 6,570) i-1, Ew(i)*aucm, Ew0(i)*aucm, edif*aucm
        write(38,570) i-1, Ew(i)*aucm, Ew0(i)*aucm, edif*aucm
      enddo
c---------------------------------------------------------------- 
c  nlin = 1, calculate vib. constants from calcu. Ev's in Ev. 
c       = 2, calculate vib. constants from readed Ev's in Ew. 
c---------------------------------------------------------------- 
        if (nlin .eq. 2) then
	        write(38,520)
          call  calVIBconst(Ew,200)
	      if (Deinp .gt. 0.0)  write(38,530)
c-
        elseif (nlin .eq. 1) then
	        write(38,540)
          call  calVIBconst(Ev,200)
	        write(38,550)
        endif
	    kw = 1
c================================================================
 500  format(//8x,"Einp(v) == INPUT energies Ev's ; ",
     #/7x,"Edif_1(v) =   Einp(v)  -  Einp(v-1) ",
     #/7x,"Edif_2(v) =  Edif_1(v) - Edif_1(v-1) ",
     #//5x,'v     Einp(v; a.u.)   Edif_1(v; a.u.)  Edif_2(v; a.u.)',/)
 510  format(//5x,'v     Einp(v; cm-1)   Edif_1(v; cm-1) ',
     #' Edif_2(v; cm-1) ',/)
 520  format(///8x,'===== Following are generated from',
     #" readed Ev's =====",)
 530  format(/6x,'===== Finish calculations for',
     #" READED Ev's =====",/)
 540  format(///6x,'===== Following are generated from',
     #" CALCULATED Ev's =====",)
 550  format(/6x,'===== Finish calculations for',
     #" CALCULATED Ev's =====",)
 570  format(3x,i3,1PE18.9,1PE18.9,x,1PE13.5,x,1PE13.5)
c----------------------------------------------------------
 900    return
      END
C===
      subroutine  guessFn(We1)
c-----------------------------------------------------------------
      implicit real*8(a-h,o-z)
C     implicit real*16(a-h,o-z)
      common /spectra/ amu,Re,De,Wee,We,WeXe,WeYe,WeZe,WeTe,WeSe,WeRe
      common /fms/ mtyp,ms,meny,Mryd,Mwe,Ny,nw0,nvd,kn,imv,kin,nlin
      common /vnumpot/ ff2, betap
      common /hnk/ hs2,Nryd,ksh
      common /fncom/ xmean(30),gf(30)
      common /fnco1/ gg1(300),gg2(300),ge1(300),ge2(300)
      common /fnco2/ gg(300),ge(300),ggf(300)
c-----------------------------------------------------------------
c  Calculate force constants ff == gg  using
c  the definitions :   f_n = [ d~nV(R) ]/[dR~n] 
c  and  the function code "dfridr" .
c
c    hs -- an estimated initial stepsize; it needs
c not be small, but rather should be an increment
c in R0 over which function fpot changes substantially.
c   err -- An estimate of the error in the derivative. 
c-----
c  Calculate numerical "force constants" for :
c    V_morse(R),  as  ks=1;     V_rydberg,  as  ks=2;
c    V_pg(R),     as  ks=3.
c
c   All force constants are the functions of De  since
c potentials V===V(R;De).
c-------------------------------------------------------- 
          ff2 = amu*We1*We1  
        betap = dsqrt(ff2/(2.0*De) )
c-------------------------------------------------------- 
           do ks=1,3
c            do k=2,Nryd+1
             do k=1,Nryd+1
               if (ks .eq. 1) then
		       gg1(k) = fmorse(k)
C=                   if (k .le. 30) gf(k) = gg1(k)
                 if (k .le. 30 .and.ksh.eq.1) gf(k) = gg1(k)
                   ge1(k) = 0.0
                     if (Ny .eq. 1) then
                       gg(k) = gg1(k)
                       ge(k) = ge1(k)
                     endif
               elseif (ks .eq. 2) then
		     ggf(k) = fryd(k,Re,We1)
                   if (Ny .ge. 2 .and. Ny .le. 25) then
                     gg(k) = ggf(k)
C=                     if (k .le. 30) gf(k) = ggf(k)
                 if (k .le. 30 .and.ksh.eq.2) gf(k) = gg(k)
                     ge(k) = 0.0
                   endif
               elseif (ks .eq. 3) then
		     gg2(k) = dfridr(k,ks,Re,hs2,err,We1)
                 ge2(k) = err
                   if (Ny .eq. 26) then
                     gg(k) = gg2(k)
                     ge(k) = ge2(k)
                   endif
                 if (k .le. 30 .and.ksh.eq.3) gf(k) = gg2(k)
               endif
             enddo
c
           enddo
c--------------------------------------------------- 
 100    return
      end 
C===
