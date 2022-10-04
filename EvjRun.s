#!/bin/tcsh
#                    EvjRUN.s
#
# A.> To evaluate energies Ev's/Evj's for diatomic molecules AB :
#  1.)  Vibrational energies Ev(i) &/or vib.-rot. energies Evj(i,j);
#  2.)  Vibrational &/or vib.-rot. THRESHOLD energies;
#  3.)  Vibrational &/or vib.-rot.  CHANNEL  energies.
#
# B.> To evaluate molecular DISSOCIATION energies De for diatomics.
#
#         !!!   Read  file  "READecm1*.fn"  for ECM1   !!!
#==========================================================================

  set   temp = '/home3/feng/temp/job1'

  cd   $temp
  rm   $temp/*

#-----
  set source = '/home1/feng/bound/AM/source'

  set target = '/home/shu/shuwork/DF/x1sigma+'
 
  set   save = '/home1/feng/bound/AM/result'

#-----
# set  mcal = 1  #  For Ev's & Evj's.
  set  mcal = 2  #  For molecular DISSOCIATION energy De.

#-----
#    As mcal = 2 , parameters 'bryd, Mryd, Ny, ak, lda' will
#           NOT be used in the calculations. You must INPUT Ev's OR
#           Delta_Ev's, and/or constant  We  in script "vibrotE.s".
#-------------------------------------------------------------------------
 
         cp  $source/vibrotE.x     .

         cp  $source/Error.fx      .
 
         cp  $target/vibrotE.s     .
 
         cp  $target/deltaG.exp    fort.23 
#-------------------------------------------------------------------------
#
#  Use following V(R) to calculate n_th force constants f_n's.
#      |
#      | = 1, use Vmorse as V(R).
#      |
#      | = 2, use  Vryd(R) = -De *[ 1 + a*x ] *exp(-a*x)         1.)
#      |                a = b * dsqrt( amu/De ) ;   x = R - Re.
#      |                b = We + bryd ;  originally, b=We.
#      |--------------------------------------------------------
#      | The GENERALIZED Rydberg potentials are defined as :
#      |
#      | = 3, use
#      |      func(R)=[exp(-a*x)/beta - (a+bryd)*x]*exp(-a*x)
#      |      Vryd(R)= De*beta*func(R)                           2.)
#      |           a = We * dsqrt( amu/De )
#      |
#      |------------------
#      | = 4, use
#      |      Uryd(R)=-De*[1.0 + (a-bryd)*x]*exp(-a*x)           3.)
#      |         a = bryd + dsqrt(bryd*bryd + f2a/De)
#      |
#      |------------------
#      | = 5, use
#      |      Uryd(R)=-De*[1.0 + a*x]*exp[-(a-bryd)*x]           4.)
#      |         a = dsqrt(bryd*bryd + f2a/De)
#      |------------------
#      | = 6, use
#      |      Umos(R)=De*[exp(-2*B*x) - (2+bryd*x)*exp(-B*x)]    5.)
#      |      Vcomb(R) = 0.5*[ Uryd(R) + Umos(R) ]
#      |             B = We * dsqrt{ amu/(2*De) }
#      |
#      |------------------
#      | = 7, use
#      |      Vmos(R)= De*[exp(-2*B*x) - 2*exp(-B*x)]
#      |      Vryd(R)=-De*[1.0 + a*x]*exp(-a*x)
#      |      V(R)=Vryd(R) + b0*[ Vryd(R) - Vmos(R) ]            6.)
#      |-----
#      | = 8, use
#      |      V(R)=Vmos(R) + b0*[ Vmos(R) - Vryd(R) ]            7.)
#      |        where  b0=bryd**Mryd  or  b0=(bryd/Re)**Mryd
#      |          and [...] has the same sign as the 1st term.
# Ny = |------------------
#      | = 9, use
#      |      Va(R)=Vcomb(R) + b0*[ Vryd(R) - Vmos(R) ]          8.)
#      |
#      |------
#      | = 10, use
#      |      Vb(R)=Vcomb(R) + b0*[ Vmos(R) - Vryd(R) ]          9.)
#      |
#      |------------------
#      | = 11, use
#      |      U(R) = [ Um(R) + Ur(R) ]/2.0                      10.)
#      |  where
#      |    Um(R) = De*[exp(-2*bt*x) - (2+b*x)*exp(-bt*x)]
#      |    Ur(R) = -De*[1 + (af-b)*x]*exp(-af*x)
#      |    bt = 0.5*(-cb + dsqrt(cb**2 + 2.0*ff2/De) )
#      |    cb = 0.5*( 1.0-sqrt(2.0) )*b ;  b == bryd
#      |
#      |------
#      | = 12, use
#      |      U(R) = [ Um(R) + Ur(R) ]/2.0                      11.)
#      |  where
#      |    Um(R) = De*[exp(-2*bt*x) - (2+b*x)*exp(-bt*x)]
#      |    Ur(R) = -De*[1 + (af-b)*x]*exp(-af*x)
#      |    bt = 0.5*(-cb + dsqrt(cb**2 + 4.0*ff2*bk/De) )/bk
#      |    cb = 0.5*( 1.0 - ak)*b  ;     b == bryd ;
#      |    bk = 1.0 + ak*ak/2.0 ;        af = ak*bt .
#      |    ak = INPUT value which is close to sqrt(2.0)=1.4142
#      |
#      |------
#      | = 13, use
#      |      U(R) = [ Um(R) + Ur(R) ]/2.0                      12.)
#      |  where
#      |    Um(R) = De*[(1-b*x)*exp(-2*bt*x) - 2*exp(-bt*x)]
#      |    Ur(R) = -De*[1 + (af-b)*x]*exp(-af*x)
#      |    bt = 0.5*(-cb + dsqrt(cb**2 + 2.0*ff2/De) )
#      |    cb = ( 1.0 - sqrt(2.0)/2 )*b  ;     b == bryd .
#      |
#      |------
#      | = 18, use Eq.12) & PVM procedures to calculate the
#      |       INItial guesses of the FORCE constants fn's
#      |       by using VIBrational constants from AM: A*X = E.
#      | 
#      |==================
#      | = 26, use Vpseudo_gaussian  as V(R).
#
# bryd - The variational constant used to adjust Vryd(R).
# Mryd = M = powers in Eq 7.) for Ny = 8.
#
#   bryd, Mryd are meaningless for Vmorse & Vpseudo_gaussian.
#
#     Code sets  beta=1.0  internally !
# beta = WIDTH (adjustable) parameter of the potential.
#
#------------------------------------------------------
#   Ny = 4, 5, 11, 12, 13, 18  are GOOD choices !
#      = 4, 5       --> f_1 =/= 0.0 which is NOT PHYsical ;
#      =11,12,13,18 --> f_1 === 0.0 which is CORRECT !
#
#   Ny = ANY values to cal.  De  as  mcal = 2 !
#
#  set Ny = 1
#  set Ny = 4    #  b = -0.379 
#  set Ny = 5    #  b =  0.379 
#  set Ny = 11   #  b = -5.13
#  set Ny = 12   #  b = -5.13 
#  set Ny = 13   #  b =  6.5252
#Nodified by Chunjun Shu
#  set Ny = 18   #  b =  6.5252
   set Ny = 18   #  b =  6.5252

#  set beta = 1.0

#-------------------------------------------------------------------------
   if ($mcal == 2) then
     echo "                                         " >! fort.43
     echo " For molecular DISSOCIATION energy De :  " >> fort.43
     echo "                                         " >> fort.43
     echo " Best De should corresponds to a best    " >> fort.43
     echo " set of VIBrational constants & a VIB.   " >> fort.43
     echo " spectrum which has a SMALL |Ave_dif|    " >> fort.43
     echo "                                         " >> fort.43
     echo "  |Ave_dif| = (SUM_v Dif_v)/(v_max+1)    " >> fort.43
     echo "                                         " >> fort.43
     echo "    Dif_v = | Einp(v) - Ecal(v) |        " >> fort.43
     echo "                                         " >> fort.43
     echo "      De = De(We,WeXe,WeYe,WeZe,...)     " >> fort.43
     echo " Ecal(v) = Ev(w0,We',WeXe,WeYe,WeZe,...) " >> fort.43
     echo "                                         " >> fort.43
     echo "    De(a.u.)           |Ave_dif|         " >> fort.43
     echo "                                         " >> fort.43
   endif
#-------------------------------------------------------------------------
#   Input the eigenvalue lda of the z component Lz of the
# electronic angular momentum L whose z axis coinciding
# with the molecular axis.
#
   set lda = 0
#  set lda = 1
#
#  mset -- Number of solution sets in A*X = E.
#                    mset >= 2   !
#  set mset = 7 
   set mset = 6 
#  set mset = 2 
   set  ikk = 50
      @ ikk += $mset
#
#  mcst -- Number of VIB. states used to solv  A*X = E.
#  set mcst = 7 
   set mcst = 8     # mcst > 7    mset + mcst =< mv
#  set mcst = 6     # mcst < 7  You'll get FEWER VIB. constants.
#-------------------------------------------------------------------------
#   bryd = ANY values to cal.  De  as  mcal = 2 !
#
# foreach  bryd (6.51 6.5252 6.5379 9.375)
# foreach  bryd (-5.13 )        #  Good for Ny=11,12
# foreach  bryd (6.5379)        #  Good for Ny=13 & De=0.173865
## Modified by Chunjun Shu
# foreach  bryd (6.5252)        #  Good for Ny=18 & De=0.17447
  foreach  bryd (6.5252)        #  Good for Ny=18 & De=0.17447

#      For  Ny=12 only ;  ak .=. dsqrt(2.0)=1.4142 :
#  foreach  ak (1.41 1.415 1.418)     
   foreach  ak (1.4142)     

#      For molecular DISSOCIATION energies De :
#               For H2 ground state, De=0.17447
#   foreach  De0 (0.170 0.1744 0.17442 0.17443 0.17445 0.176 0.1863)
#   foreach  De0 (0.1696865)
#   foreach  De0 (0.2250829)
    foreach  De0 (49383.53)
#   foreach  De0 (0.17447)
#   foreach  De0 (0.1744286)
#   foreach  De0 (0.1737238)
#   foreach  De0 (0.0)

      nice +4  vibrotE.s   $bryd $De0 $Ny $ak $lda $mcal $mset $mcst

      cp  fort.6           Out.De$De0
#     cp  fort.6     $save/Out.De$De0
#     cp  fort.6     $target/data/Out.Evj$bryd

      cp  fort.31       Wxyztsr.B$bryd
#     cat fort.31 >> $target/Wxyztsr.B$bryd

      cp  fort.33      Fn1to8.B$bryd
#     cat fort.33 >> $target/Fn1to8.B$bryd

      cp  fort.35       ConstEvj.B$bryd
#     cp  fort.35       ConstEvj.ak$ak.B$bryd
#     cp  fort.35    $target/data/ConstEvj.B$bryd

      cp  fort.38       CompEvDe$De0
      cp  fort.38       EvDe$De0.B$bryd
#     cp  fort.38    $target/data/CompEvDe$De0.B$bryd

      cat fort.40  >>   fort.43
      cat fort.41  >>   fort.44
      echo "    "  >>   fort.44

#   ik1 -- Starting unit to write information for De. 
#   ik2 -- Starting unit to  SAVE information for De. 

        set ik1 = 51
        set ik2 = 71

      doloop:
        if ($ik1 <= $ikk) then
          cat fort.$ik1 >>   fort.$ik2
             @ ik1++
             @ ik2++
          goto  doloop
        endif
#
        if ($mcal == 1) then
          cp    fort.33    fn.B$bryd
        endif

    end
   end
  end
#
           set ik1 = 51
           set ik2 = 71

        echo "    "  >>   fort.44
        echo "       LOOPS over different set : " >> fort.44
        echo "    "  >>   fort.44

      doloop1:
        if ($ik1 <= $ikk) then
	      cat fort.$ik2  >>   fort.44
              echo "    "  >>   fort.44
             @ ik1++
             @ ik2++
          goto  doloop1
        endif 
#
	    cp  fort.43             DeCMP1.B$bryd
	    cp  fort.44             DeCMP2.B$bryd
#         cat fort.44  >> $target/DeCMP2.B$bryd
#
#==========================================================================

