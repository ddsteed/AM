C    program Err
C========================================================================
      IMPLICIT  REAL*8(A-H,O-Z)
      
      Err = 5.000000E-4
      Deper = -1.000E-4
      open(unit = 5, file="Ave.error.tmp")
      open(unit = 8, file="Dediff.tmp")
      open(unit = 6, file="Ave.error.Ok")
      open(unit = 7, file="Dediff.Ok")

      Do 110 I = 1,6
      read(8,*) Dediff
      read(5,*) AveErr

      if((Dediff.LE.-0.00001).AND.(Dediff.GE.Deper)) then
C      read(5,*) AveErr

      if(AveErr.LE.Err) then
      write(6,*) "T"
      write(7,*) "T"
      goto 10
      endif

      endif

110   CONTINUE
      close(5)
      close(8)
      close(6)
      close(7)
10    CONTINUE
      end
