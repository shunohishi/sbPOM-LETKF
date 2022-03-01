          program windstress

          implicit none
          real wvel, wind, cd
          integer i

          do i=0,1000
      
          wvel=i*0.1

C Large and Pond (1981) 
          if(wvel.lt.10.) then
            cd=1.14
          else if( wvel.lt.26. ) then
            cd=0.49+0.065*wvel
          else
            cd=0.49+0.065*26.
          end if

C Garrat (1977)
          cd=0.75+0.067*wvel
          if( wvel.gt.26. ) then
            cd=0.75+0.067*26.
          end if

C Mellor and Blumberg (2004)
          if( wvel.lt.3.7313 ) then
            cd=1.5-wvel*0.134002
          else
            cd=0.75+0.067*wvel
            cd=min(2.492,cd)
          end if

          write(6,*) wvel, cd

          end do

          end
