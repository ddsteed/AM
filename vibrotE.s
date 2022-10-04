#!/bin/tcsh
#                 vibrotE.s
#-----------------------------------------------------------------------
# A.> To evaluate ENERGIES Ev's/Evj's for diatomic molecules :
#  1.)  Vibrational energies Ev(i) &/or vib.-rot. energies Evj(i,j);
#  2.)  Vibrational &/or vib.-rot. THRESHOLD energies; 
#  3.)  Vibrational &/or vib.-rot.  CHANNEL  energies. 
#
# B.> To evaluate molecular DISSOCIATION energies De for diatomics :
#                 set  mtp = 3
#-----------------------------------------------------------------------

@ mlsta = 6
set   temp = '/home/shu/shuwork/DF/x1sigma+/result'
set   target = '/home/shu/shuwork/DF/x1sigma+'
if(-s $temp/data) then
rm $temp/data
endif

nextshu:


# set mtp = 1   # use input vib-rot constants to cal. Ev & Evj
#                   Set nlin=0 for mtp=1 !
# set mtp = 2   # use f_n's -- derivatives of V(R) to cal. Ev & Evj
  set mtp = 3   # calc. DISSOCIATION energy De from input {Ev}.

  set mvr = 0   # for both vibration and vib.-rotation
# set mvr = 1   # for vibration 
# set mvr = 2   # for vib.-rotation

  set  ne = 1   # No. of scattering energies for which the scattering 
                # threshold and channel energies are wanted.

  set  nv = 27  # nv = number of vibrational states considered

  set rmx = 13.0  #  Checking for  V(R->rmx=inf) = De .
  set nry = 80    #  No of f_n's used to cal. V(R->rmx=inf). 

  set bryd = $argv[1]
  set  De0 = $argv[2]
  set   Ny = $argv[3]
  set   ak = $argv[4]
  set  lda = $argv[5]
  set mcal = $argv[6]
  set mset = $argv[7]
  set mcst = $argv[8]

 if ($mcal == 2) then
  set  mtp = 3   # use INPUT Ev's OR Delta_Ev's to calculate De.
 endif

#  mtd -- Decide which VIB. states being used to solve A*X = E :
#       = 0, use VIB. states in successively increasing order;
#       > 0, use VIB. states particularly specified.
#
# set  mtd = 0  
# set  mtd = 1   # for nw0 = 0;  mv = 14;  mcst = 7;  mset = 7.  
  set  mtd = 2   # for nw0 = 1;  mv = 14;  mcst = 8;  mset = 6.
#=======================================================================
#   All reading are FREE formats !   For H_2 X_1sigma_g+ state.
#-----------------------------------------------------------------------
#         mtyp
   echo " $mtp "  >! fort.5   
#---
#  nw0 = 2, E(v) = w0 + (We + we0)*(v+0.5) - WeXe*(v+0.5)**2 + ...
#               We from INPUT;  w0 & we0 from perturbation cal.
#      = 1, E(v) = w0 + We'*(v+0.5) - WeXe*(v+0.5)**2 + ...
#               We' from solving  A*X = E .
#      = 0, E(v) = We*(v+0.5) - WeXe*(v+0.5)**2 + ...
#  nvd = number of vib. states used in comparison with
#          the input Ev or energy differencies. nvd =< n_expt.
#          nvd is slightly smaller than n_expt which is the
#          no. of input Ev OR energy differencies for comparison.
#             The mv below must be  mv =< nvd ;    nvd =< 95.
#   kn = number of EXTRA vib. states used in comparison.
#  imv = maximum vib. quantum number used to calculate De(inp)
#          using INPUT vibrational constants.
#        imv is used in subroutine FindnumDe only.
#  kin = 1, READ expt'l energies or their differencies; 
#      = 0, do NOT read.
# nlin = 1, cal. fn's from CALC. vib. energies using "broydn";
#      = 2, cal. fn's from INPUT vib. energies using "broydn";
#      = 0, do NOT cal. force constants fn's.
#   Set  nlin = 0  for mtp = 1 !  Otherwise it may cause ERROR !
#
#          nw0 nvd  kn   imv kin nlin
#  echo  "  2   13  20  1900  0   1  "  >> fort.5
   echo  "  1   $nv 15  1900  1   2  "  >> fort.5
#  echo  "  0   15  20  1900  1   2  "  >> fort.5
#  echo  "  0   11  20  1900  1   2  "  >> fort.5
#  echo  "  1   14  20  1900  1   0  "  >> fort.5
#---
# mset -- number of set of solution vectors  X;
#           usually, set  mset =< 20.
#           !  mset =<  mv - mcst  !   for mtd = 0
#
# mcst -- number of vib. constants in a vector X;
#           usually, set  mcst =< 7.  BUT as solve
#            for X using INPUT Ev's,  7 =< mcst < mv.
#               A(mcst, mcst);  b(mcst, 1)
#
#   mv -- number of vib. states in A used to
#           solve linear equation    A * X = b
#             set  mv =< nvd ;    nvd =< 95.
#    X = ( We, WeXe, WeYe, WeZe, WeTe, WeSe, WeRe, ...)
#
# meny = 0 :
#  Ev = We(v+1/2) - WeXe(v+1/2)**2 + WeYe(v+1/2)**3
#                 + WeZe(v+1/2)**4 + WeTe(v+1/2)**5 + ... + ...
# meny = 1 :
#  Ev = We(v+1/2) - WeXe(v+1/2)**2 + WeYe(v+1/2)**3
#                 - WeZe(v+1/2)**4 + WeTe(v+1/2)**5 - ... + ...
#  
# mdp -- The mdp_th set from (A*X=E) CONVERGES spectrum { Ev }.
#  
#          mset  mcst  mv meny  mtd mdp
  echo  " $mset $mcst  $nv 0   $mtd  5 "  >> fort.5
# echo  "    4     7   14  1   $mtd  3 "  >> fort.5
#---
#         For  mtd = 0 :
# nvs(i,1) -- value of INITIAL vib. quantum state
#             v+1 in the ith set of mset vectors, e.g.
#               1  = 1, 2, 3, 4, 5   (mcst=5 states)
#         nvs(i,1) = 3, 2, 6, 8, 9
#            v_ini = 2, 1, 5, 7, 8
#        ! nvs(i,1) + mcst =< mv   !   
#
#    When mcst=6 & mtd=0 , the states used in a set are :
#  nvs(1,1)=3 :  v = 2, 3, 4, 5, 6, 7.
#  nvs(3,1)=6 :  v = 5, 6, 7, 8, 9,10.
#
#    For ith set, if mcst=7 & mtd > 0 , ARBItrary states are used :
#  nvs(i,j)=2, 3, 6, 7, 14, 15.
#         v=1, 2, 5, 6, 13, 14.
#
#              Specify INITIAL VIB. state !
 if ($mtd == 0) then
#     For mtd = 0 :  nvs(i,1), i=1,mset
# echo  "  1   2   3   4   5   6   7   8 "  >> fort.5
# echo  "  1   2   3   4   5   6   7 "  >> fort.5
# echo  "  1   2   4   6   8  10  "  >> fort.5
##echo  "  1   2   4   5   6   7  "  >> fort.5
  echo  "  3   2   4   5   6   7  "  >> fort.5
# echo  "  3   2   4   6   5   7  "  >> fort.5
# echo  "  3   2   4   6          "  >> fort.5
# echo  "  1   2                  "  >> fort.5
 endif
#----          Specify EVERY VIB. state !
#                        nvs(i,mcst) =< mv  !
 if ($mtd == 1) then
#     For mtd = 1 :  nvs(i,j): i=mset; j=1, mcst
# echo  "  1   2   3   4   5   6   7   8   9  10 "  >> fort.5
# echo  "  1   2   4   6   8  10  13  15  "  >> fort.5 
# echo  "  1   2   4   7   10 13  14  15  "  >> fort.5 
#-
  echo  "  2   4   6   8  10  13  15  "  >> fort.5  # BEST !
  echo  "  2   3   6   8  10  13  15  "  >> fort.5  
  echo  "  2   3   5   7  10  13  15  "  >> fort.5
  echo  "  2   3   5   7   9  11  12  "  >> fort.5
  echo  "  1   3   6   8  10  13  15  "  >> fort.5  
  echo  "  1   2   5   7  10  13  15  "  >> fort.5
  echo  "  1   2   5   8  10  13  15  "  >> fort.5  
# echo  "  1   2   4   7   9  12  14  "  >> fort.5 
# echo  "  1   3   5   7   9  11  13  "  >> fort.5
 endif
#---
 if ($mtd == 2) then
#     For mtd = 2 :  nvs(i,j): i=mset; j=1, mcst
head -$mlsta $target/msetresult.dat >! tempshu

 set AveE = `tail -1 tempshu`
 if ("$AveE" == "End") then
 tail -7 tempshu  >! tempshu1
 head -6 tempshu1 >> fort.5
 else
 tail -6 tempshu >> fort.5
 endif


# echo  "  2   4   6   8  10  12  13  15  "  >> fort.5  # 
# echo  "  2   3   6   8  10  11  13  15  "  >> fort.5  
# echo  "  2   3   5   7  10  11  13  15  "  >> fort.5
# echo  "  2   3   5   7   9  11  13  14  "  >> fort.5
# echo  "  1   3   6   8  10  12  14  15  "  >> fort.5  
# echo  "  1   2   5   7  10  12  14  15  "  >> fort.5
 endif
#---
#   The initial guess for vibrational FORCE constants :
#    (They may be used to cal. gf(i) using 'broydn')
#      gf(i) < 0, i = 3, 5, ...;  gf(i) > 0, i = 4, 6, ...
#              gf(i): i = 3, mset+2
  echo  " -0.8  1.2 -1.5  1.7 -1.9  2.0 -2.1 2.3 "  >> fort.5
#==
#    When check VIBrational energies using the vibrational
#  constants generated from the INPUT vib. energies by 
#  linear algebraic method (A*X=B).
#  bxe, ..., bre = 0, Zero calc. constants :
#                    WeXi, WeYi, WeZi, WeTi, WeSi, WeRi;
#                = +1, Do NOT change sign of these constants;
#                = -1, Change the sign of these constants.
#         bxe  bye  bze  bte  bse  bre
  echo " -1.0  1.0  1.0  1.0  1.0  1.0  "  >> fort.5
# echo "  1.0  1.0  1.0  1.0  1.0  1.0  "  >> fort.5
#---
# Mvr  = 1, for vibration; = 2 for rot-vibration; =0 for BOTH.
#  NE  = The # of input SCATTERING energies
#  Nv  = # of vibrational states considered
#            Nv=6 means that v=0,1,2,3,4,5
# iniv = Quantum # of the initial vibrational referen#e state
# jnir = Quantum # of the initial  rotational referen#e state
#                 (Usually take iniv=v0=0, jnir=j0=0)
#  Nw  = 0, input (WE,...,Be,...) in a.u.;
#      = 1, input (WE,...,Be,...) in cm-1, code will convert them.
#
#          Mvr   NE   Nv  iniv jnir  Nw
#  echo " $mvr  $ne  $nv    0    0    0  "  >> fort.5   
   echo " $mvr  $ne  $nv    0    0    1  "  >> fort.5   
#  echo " $mvr  $ne  $nv    0    0    1  "  >> fort.5   
#---
#   Re = The equilibrium internuclear distance.
#   De = The dissociation energy for AB.
#          If Nw=0, Re in ao, De in a.u.; or Re in A, De in cm-1. 
#
#      INPUT De = 0.0, if you were to cal. the UNKNOWN De.
#
#   Code sets  beta = 1.0  internally !
# beta = WIDTH (adjustable) parameter of the potential.
#           Re      De       
#  echo " 1.7327238 $De0     " >> fort.5  
   echo " 0.916930  $De0     " >> fort.5  
#  echo " 1.4011  0.0       " >> fort.5  
#  echo " 1.4011  0.174457  " >> fort.5  # De in a.u.
#  echo " 1.4011  0.173865  " >> fort.5  # De <- Kolos (1968)
#  echo " 1.4011  38289.54  " >> fort.5  # De in cm-1
#---
# ams = The mass for atom A in atomic_mass_unit (amu).
# bms = The mass for atom B in amu.
# am0 = The reduced mass for molecule AB in amu.
#   If input am0 = 0.0, code will calculate am0.
#           ams        bms        am0
#  echo " 1.0078252  6.0151200    0.0   " >> fort.5 
   echo " 2.0141018  18.998400    0.0   " >> fort.5 
#---
#   nj(i) = # of rotational states in the i_th vibrational state.
#     (nj-1) - the HIGHEST rotational state in v state.
#       To have  [ SUM of Evj(i,j) ] = E(v_max) = De :
#     Set nj >= 1 + j_max (of vib. state) = 1 + (nj-1)
#               nj(i) = 1, nv
#  echo  "  9  9  9  9  9  9  9  9  9  9 "  >> fort.5 # mtp=1
   echo  "  9  9  9  9  8  8  7  6  5  5 "  >> fort.5 # mtp=2
 if ($nv > 10) then
#  echo  "  9  8  8  8  9  8  7  6  6  5 "  >> fort.5 # mtp=1
   echo  "  4  3  2  1  1  1  1  1  1  1 "  >> fort.5 # mtp=2
 endif
 if ($nv > 20) then
   echo  "  4  2  2  1  1  1  1  1  1  1 "  >> fort.5 # mtp=1
 endif
 if ($nv > 30) then
   echo  "  4  2  2  1  1  1  1  1  1  1 "  >> fort.5 # mtp=1
 endif
#---
#   Input the eigenvalue lda of the z component Lz of the
# electronic angular momentum L whose z axis coinciding
# with the molecular axis.  "lda" presents the coupling
# between the electronic and nuclear motion of the diatomic
# system.
#          LAMDA
   echo "  $lda  " >> fort.5  
#---
#            E(i; eV)   i=1,NE
   echo " 1.0  2.0  3.0  4.0  5.0  " >> fort.5  
#
#--------------------------
#       Vibrational & rotational constants :
#   When Nw=0, input them in a.u., or in cm-1.
#              We           WeXe          WeYe      
#  echo "  2.0062460D-02  5.528202D-04  3.703845D-06  " >> fort.5 
#-     G. Herzberg (1979) :
#  echo "  6.4685600D-03  5.528475D-04  3.703845D-06  " >> fort.5 
   echo "  2.9981920D+03  4.576100D+01  0.0           " >> fort.5 
#-
#            WeZe         WeTe          WeSe 
#  echo " -2.278167D-08    0.0           0.0   " >> fort.5 
   echo "  0.0             0.0           0.0   " >> fort.5 
#-
#            WeRe        
#  echo "    0.0     " >> fort.5 
   echo "    0.0     " >> fort.5 
#---
#             Be          Alpha_e        gammae
#  echo "  2.772666D-04  1.395150D-05   2.597111D-07 " >> fort.5
   echo "  1.101020D+01  3.017000D-01   0.0          " >> fort.5
#-
#            Der          betae
#  echo "    0.0           0.0   " >> fort.5
   echo "    0.0           0.0   " >> fort.5
#------------------
#    When calculate VIB-ROTational energies for mtyp=2 :
#  aye, ..., are = 0, Zero calc. constants 
#                    WeY, WeZ, WeT, WeS, WeR;
#                = +1, Do NOT change sign of these constants;
#                = -1, Change the sign of these constants.
#            aye  aze  ate  ase  are
     echo "  1.0  1.0  1.0  1.0  1.0  "  >> fort.5
#---
#  abe,aae,age,
#  ae3 --- ae7 = 0, Zero calc. ROT. constants Eta_3 -- Eta_7;
#              = +1, Do NOT change sign of these constants;
#              = -1, Change the sign of the cal. ROT. constants.
#            abe  aae  age  ae3  ae4  ae5  ae6  ae7 
     echo "  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  " >> fort.5
#      abe -> Be;  aae -> Alpha_e;  age -> Gamma_e.
#---
#  ade,abt,ax2,
#  ax3 --- ax7 = 0, Zero calc. ROT. constants Xsi_2 -- Xsi_7;
#              = +1, Do NOT change sign of these constants;
#              = -1, Change the sign of the cal. ROT. constants.
#            ade  abt  ax2  ax3  ax4  ax5  ax6  ax7 
     echo "  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  "  >> fort.5
#    echo " -1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  "  >> fort.5
#
#      ade -> D_e;  abt -> Beta_e.
#      ade =  1.0, -> + D_e;  ade = -1.0, -> - D_e .
#    A. L. G. Rees, Proc. Phys. Soc.(London)59, 998(1947) :
#              E(v,j) = ..., + D_e*[J*(J+1)]**2 + ...
#    G. Herzberg (1953) :
#              E(v,j) = ..., - D_e*[J*(J+1)]**2 + ...
#--------------------------
 if ($mtp >= 2) then  
#
#  Use following V(R) to calculate n_th force constants f_n's.
#      |
#      | = 1, use Vmorse as V(R).
#      |
#      | = 2, use  Vryd(R) = -De *[ 1 + a*x ] *exp(-a*x)         1.)
#      |                a = b * dsqrt( amu/De ) ;   x = R - Re.
#      |                b = We + bryd ;  originally, b=We.
#      |
#      |====================================================
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
#      |-----
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
#      |
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
#      |    cb = 0.5*( 1.0-sqrt(2.0) )*b ;  b == bryd .
#      |
#      |------
#      | = 12, use
#      |      U(R) = [ Um(R) + Ur(R) ]/2.0                      11.)
#      |  where
#      |    Um(R) = De*[exp(-2*bt*x) - (2+b*x)*exp(-bt*x)]
#      |    Ur(R) = -De*[1 + (af-b)*x]*exp(-af*x)
#      |    bt = 0.5*(-cb + dsqrt(cb**2 + 4.0*ff2*bk/De) )/bk
#      |    cb = 0.5*( 1.0 - ak)*b  ;     b == bryd  ;  
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
#      |=====================================
#      | = 26, use Vpseudo_gaussian  as V(R).
#
#    Ny = 4, 5, 11, 12, 13, 18  are GOOD choices !
#       = 4, 5       --> f_1 =/= 0.0 which is NOT PHYsical ;
#       =11,12,13,18 --> f_1 === 0.0 which is CORRECT !
#
#   Ner = number of force constants used to cal. ro-vib constants.
#
#  bryd - The variational constant used to adjust Vryd(R).
#  Mryd = M = powers in Eq 7.) for Ny = 8 .
#
#      bryd, Mryd are meaningless for Vmorse & Vp_g.
#
#  hs1, hs2 -- Estimated initial stepsize used by code "dfridr".
#                Good range for hs1 :  0.001 --> 4.0  for H_2.
#                  hs1 for Vrydberg;   hs2 for V_p-g .
#
#            Ner  Ny   bryd  Mryd   hs1   hs2   
     echo  "  7  $Ny  $bryd    1   2.000  2.50  " >> fort.5
#    echo  "  7  $Ny  $bryd    2   2.000  2.50  " >> fort.5
#---
#    V(R->inf=Rmax) = De = Sum_n  f_n*(Rbig - Re)**n/n!
#   Mwe = 1, New  We = We + bryd ;  = 0, We = We. 
#  Nryd = Number of f_n's cal. using Vryd(R).
#  Rmax = The maximum R value for  V(R->Rmax=inf).
#   ak = a value close to sqrt(2.0)=1.4142 for Ny=21.
#         ak is used in alpha = ak*beta  in potential :
#           U(R) = [ Um(R) + Ur(R) ]/2.0
# Rless = The R value in    Rbig = Rmax +/- Rless.
# Dconv = The tolerance used to check if V(Rmax) = De.
#
#             Mwe  Nryd  Rmax   ak  Rless  Dconv
     echo  "   0   $nry  $rmx  $ak  0.01   5E-08  "  >> fort.5
#    echo  "   1   $nry  $rmx  $ak  0.01   5E-07  "  >> fort.5
#---
# br4,br6,br8,br10,br12,br14 = The scaling switches for the
# long range force const f_n(long;R=Re) = d~n/(dR~n)[1/R**k]
# which will be added into f_n(R=Re) = V_rydberg~n(R=Re).
#   f_n(long;R=Re) = (-1.0)**n *(k+n-1)!/[(k-1)!*R**(k+n)]
# which are important for Wande-Waals molecules and 
# quasi-stable molecules.
#    brn = 1.0, Add f_n(long;R=Re); = 0.0, do NOT add.
#
#             br4  br6  br8  br10 br12 br14
#    echo  "  1.0  1.0  1.0  1.0  1.0  1.0  "  >> fort.5
#    echo  "  1.0  0.0  0.0  0.0  0.0  0.0  "  >> fort.5
     echo  "  0.0  0.0  0.0  0.0  0.0  0.0  "  >> fort.5
#---
# ksh = 0, use the INPUT first guess for FORCE constants fn's ;
#     = 1, use fn's from Morse V(R);
#     = 2, use fn's from Rydberg-like V(R);
#     = 3, use fn's from P-G V(R).
#            ksh 
   if ($mtp == 2) then
     echo  "  2  "  >> fort.5
   endif
   if ($mtp == 3) then
#    echo  "  0  "  >> fort.5
     echo  "  2  "  >> fort.5
   endif
 endif 
#=======================================================================

 echo "    "
 echo "    "
#echo " Calculate De OR Ev(i) &/or Evj(i,j) for AB molecule ! "
 echo "    "

  
  nice  +4  vibrotE.x  < fort.5 >! fort.6


#echo " *** Finish vibrotE for De = $De0 &  bryd = $bryd ***  "
#echo " *** Finish vibrotE for ak = $ak  &  bryd = $bryd ***  "

# rm *
if(-s $temp/Dediff.Ok) then
rm $temp/Dediff.Ok
endif
if(-s $temp/Ave.error.Ok) then
rm $temp/Ave.error.Ok
endif
 nice +4  Error.fx
echo "Set: $mlsta"


set Desdiff  = `tail -1 $temp/Dediff.Ok`
set Dess  = `tail -1 $temp/Dediff.tmp`
echo "Dissociation energy%:\t$Dess"

if ($Desdiff == "T" ) goto Error
if ("$AveE" == "End") goto success
@ mlsta = $mlsta + 6
goto nextshu

Error:
set Ave_error = `tail -1 $temp/Ave.error.Ok`
set Ave_error_all = `tail -1 $temp/Ave.error.tmp`
echo "Average Error%:\t\t$Ave_error_all"
echo "   "

if ("$AveE" == "End") goto success


if ($Ave_error == "T" ) then

echo "$mlsta"   >> $temp/data
echo "-------"  >> $temp/data

#  set AveE = `tail -1 tempshu`
#  if ("$AveE" == "End") then
#  tail -7 tempshu  >! tempshu1
#  head -6 tempshu1 >> fort.5
#  else
 tail -6 tempshu >> $temp/okdata 
# endif

cat DeEr.all    >> $temp/data

# Modified by Chunjun Shu
#   goto success

endif

@ mlsta = $mlsta + 6
goto nextshu
success:
echo "End"  >> $temp/okdata
