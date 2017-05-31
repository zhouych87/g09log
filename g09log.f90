program g09log
  implicit none
  character(20):: fname
  character(27)::tmp1,mem,ncpu,chk
  character(79)::tmp,setting
  character(9)::val
  character(1)::arg2
  character(2)::elemnt(1000)
  integer::funct, nmo,lngth,fnmo,nline,i,j,nm
  real:: hm,lm,crd(3,1000)
  logical:: omo,debug ! 0 false occupation 
  
  debug=.false. !
!  debug=.true.
  if (iargc()/=2) then 
     write(*,*) "input parameter missing,stop"!    ch1(ii)=c(ii,imol1)
     write(*,*) "please input the name of the gaussian output file, functions. "
     write(*,*) " FOR example: pcbm.log  1 " ! cl1(ii)=c(ii,imol1+1) !! add by Yecheng Zhou
     stop
  end if 
  
  call getarg(1,fname)
  call getarg(2,arg2)
  read(arg2,*) funct ! 0 convg; 1 homo,lumo;3 log2gjf 

if (funct==0) then 
!   setting="grep 'Predicted change in Energy' "//fname
!   fname2="|sed 's/\=/\ \ /'| sed 's/D/E/'"
!   write(*,'(a,a)') trim(setting),fname2
   setting="grep 'Predicted change in Energy' "//trim(fname)//"|sed 's/\=/\ \ /'| sed 's/D/E/'"
   write(*,*) trim(setting)
  CALL SYSTEM(trim(setting))
  stop 
end if 

 open(unit=99,file=fname,status='old',err=9993)

if (funct==1) then ! obtain homo and lumo and its orbital number
  nmo=0 
 omo=.True.
 nline=0

  do while(.true.)
    read(99,'(a27)',end=9994) tmp1
    nline=nline+1
    if(tmp1 .eq.' Alpha  occ. eigenvalues --') then 
          omo=.True.
          nmo=nmo+5
        
    elseif((omo .eqv. .True.) .and. (tmp1 .eq.' Alpha virt. eigenvalues --')) then 
       backspace(99) 
       nmo=nmo-5 
       backspace(99) 
       read(99,'(a79)') tmp
       lngth=len(trim(tmp))-26 
       if (debug .eqv. .true. ) then 
              write(*,*) tmp
              write(*,*) len(trim(tmp)),len(tmp),lngth,int(lngth/10)
              write(*,*) "number of orbitals:", nmo  
        end if 
       nmo=nmo+int(lngth/10)
       val=tmp((len(trim(tmp))-8):(len(trim(tmp))+2))
       read(val,*) hm
       
       read(99,'(a79)') tmp
        if (debug .eqv. .true. ) write(*,*) tmp 
       val=tmp(30:39)
       read(val,*) lm
       omo=.false.
       fnmo=nmo
        nmo=0 
        if (debug .eqv. .true. ) write(*,*) nline 
    else if(tmp1 .eq. ' Normal termination of Gaus') then 
        exit
    end if 
 enddo

  if (debug .eqv. .true. ) write(*,"(4(f10.6x))") hm,lm,27.2*hm,27.2*lm 
  if (debug .eqv. .true. ) write(*,*) nline 
  write(*,"(2(i5x,f10.6x))") fnmo,27.2*hm, fnmo+1, 27.2*lm 
  stop

else if (funct==2) then ! funct=2, write gjf from log 
   omo =.false.
  do while(.true.)
    read(99,'(a)',end=9994) tmp
    !write(*,*) trim(tmp) 
    
    if (omo .eqv. .false.) then 
        if (tmp(1:5) .eq. ' %mem') then 
          mem=tmp(2:20) 
        else if (tmp(1:6) .eq. ' %npro') then 
          ncpu= tmp(2:20) 
        else if (tmp(1:2) .eq. ' #') then 
          setting= tmp(2:len(tmp))
          chk=fname(1:(len(trim(fname))-4))
          write(*,*) 'chk',chk,fname(1:(len(trim(fname))-4)),fname
          write(*,*) "setting ok"
        else if (tmp(1:19) .eq. ' Symbolic Z-matrix:') then 
              write(*,*) "going to read element and coordination"
            read(99,'(a)',end=9994) tmp
             i=0
             do while(.true.)
                i=i+1 
                read(99,*,err=9991, end=9995) elemnt(i),(crd(j,i),j=1,3)
              end do 
         
9991         backspace(99)
             read(99,'(a)',end=9994) tmp
             if (tmp(1:1) .eq. ' ') then 
                  write(*,*) 'Setting and element read finished'
                 nm=i-1
             else 
                 write(*,*) tmp
                 write(*,*) 'error in reading element'
             end if
            omo=.true.
             write(*,*) "setting all ok"
             write(*,*) mem,chk,ncpu  
        end if ! crd 
        
     else if (tmp(1:45) .eq. & 
     & '                          Input orientation:'.or. &
     & tmp(1:47) .eq. & 
     & '                         Z-Matrix orientation:') then 
          write(*,*) 'Read coordination once again'
       !   write(*,*) tmp
       !   write(*,*) trim(tmp) 
          
          do i=1,4
             read(99,'(a)')
          end do 
          do i=1,nm
             read(99,*,end=9996) j,j,j,(crd(j,i),j=1,3)
          end do 
     end if 
     
    if(tmp(1:27) .eq. ' Normal termination of Gaus') then 
       Write(*,*) 'Gaussian terminated normally'
        exit 
     end if 
   end do 
   
   write(*,'(a,a,a)') 'Going to write ', trim(chk), 'mo.gjf' 
   call writegjf(mem,ncpu,chk,setting,elemnt,crd,nm)
   write(*,*) 'writing finished'
   stop 
end if 
    
9993 write(*,*) "not exist"
9994 write(*,*) "Calculation is not finished or not succeed"
9995 write(*,*) "element"
9996 write(*,*) "crd"
write(*,*) tmp

   write(*,'(a,a,a)') 'Going to write ', trim(chk), 'mo.gjf'
   call writegjf(mem,ncpu,chk,setting,elemnt,crd,nm)
   write(*,*) 'writing finished'
   stop 

end program g09log 



      SUBROUTINE writegjf(mem,ncpu,chk,setting,elemnt,crd,nm)
      implicit none
      real::crd(3,1000)
      character(2)::elemnt(1000)
      character(30):: fname
      character(27)::mem,ncpu,chk
      character(79)::setting  
      integer::i,j,nm
      
      
      write(fname,'(a,a)') trim(chk), 'mo.gjf' 
      open(320, file=trim(fname),status = 'replace')
      write(320,'(a,i0)') trim(ncpu)
      write(320,'(a)')  trim(mem)
      write(320,'(a,a,a)') '%chk=',trim(chk),'.chk'
      write(320,'(a)') trim(setting)
      write(320,'(a)') ' '
      write(320,'(a)') 'MO generated'
      write(320,'(a)')  ' '
      write(320,'(a)') '0 1' ! should check it
      
      do i=1,nm
        write(320,'(a,3(3xf12.6))') elemnt(i),(crd(j,i),j=1,3) 
      end do
      write(320,*) ' '
      write(320,*) ' '
      close(320)
      END SUBROUTINE 
