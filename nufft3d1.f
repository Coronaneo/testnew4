	subroutine nufft3dI(nj,x,iflag,ns,rt,tol,U,V,xxsub,r)
	implicit none

	integer :: ns,rt,iflag,nj,k(ns*ns*ns,3),j,r,xsub(nj,2),i,xxsub(nj)

	real*8 pi,x(nj,3),tol
	parameter (pi=3.141592653589793238462643383279502884197d0)
	complex*16 fftconst,U(ns*ns*ns,r),V(nj,r)
	complex*16 M(ns*ns*ns,nj)
        integer engOpen,engGetMatrix,mxCreateFull,mxGetPr
        integer ep,T,D
        integer engPutMatrix,engEvalString,engClose
        integer eps,rm,UU,VV
        
        do w = 1,ns
	   do j = 1,ns
	      do i = 1,ns
	         k((w-1)*ns*ns+(j-1)*ns+i,1) = i
	         k((w-1)*ns*ns+(j-1)*ns+i,2) = j
                 k((w-1)*ns*ns+(j-1)*ns+i,2) = w
	      enddo
	   enddo
        enddo

	fftconst = iflag*dcmplx(0,1)/ns*2*pi

	do i = 1,ns
	   do j = 1,nj
	      M(i,j) = exp(fftconst*(k(i,1)*(x(j,1)-floor(x(j,1)+0.5))+
     &     (k(i,2)*(x(j,2)-floor(x(j,2)+0.5))+(k(i,2)*(x(j,3)-floor(x(j,3)+0.5))))
	   enddo
	enddo
   
        ep = engOpen('Matlab')
        T = mxCreateFull(ns*ns*ns,nj,0)
        call mxSetName(T,'T')
        call mxCopyComplex16ToPtr(M,mxGetPr(T),ns*ns*ns*nj)
        call engPutMatrix(ep,T)
        eps = mxCreateFull(1,0,0)
        call mxSetName(eps,'eps')
        call mxCopyReal8ToPtr(tol,mxGetPr(eps),1)
        call engPutMatrix(ep,eps)
        rm = mxCreateFull(1,0,0)
        call mxSetName(rm,'rm')
        call mxCopyIntegerToPtr(rt,mxGetPr(rm),1)
        call engPutMatrix(ep,rm)
        call engEvalString(ep,'[U,V]=lowrankfac(T,eps,rm,rm)')
        UU = engGetMatrix(ep,'U')
        VV = engGetMatrix(ep,'V')
        call mxCopyPtrToComplex16(mxGetpr(UU),U,ns*ns*ns*r)
        call mxCopyPtrToComplex16(mxGetpr(VV),V,nj*r)
	call lowrankfac(M,tol,rt,rt,U,V)

	xsub = mod(floor(x+0.5),ns)+1
	do i = 1,nj
	   xxsub(i) = xsub(i,2)*nj-nj+xsub(i,1)
	enddo

	r = size(V,2)

	end subroutine


	subroutine nufft3dIapp(nj,plan,c,U,V,xxsub,ns,iflag,r,S)
	implicit none
	integer  r,i,j,k,nj,ns,iflag,num
	integer mm
	integer xxsub(nj)
	complex*16 M(r,ns*ns*ns),N(ns,ns,ns,r),S(ns*ns*ns),c(nj),U(ns*ns*ns,r)
        complex*16 V(r,nj)
        complex*16 NN1(ns*ns,r),SS(ns,ns,ns),SS1(ns/2,ns/2),NN2(ns,ns,r)
	complex*16,allocatable :: NN(:,:,:,:)
	double complex in1, out1
	real*16  time1,arr(4)
	dimension in1(ns,ns,ns), out1(ns,ns,ns)
	integer*16 :: plan
        character*8 date
        character*10 time
        character*5 zone 
        integer*4 values1(8),values2(8)

        arr(1)=3600
        arr(2)=60
        arr(3)=1
        arr(4)=0.001
        num=100
        !print *,'ok1',size(xxsub)
        !print *,'xxsub=',size(xxsub),xxsub
        call date_and_time(date,time,zone,values1)
        do j = 1,num
ccccccccccccccccccccccccccccccccccccccccccccccccccc
        !part 1
	M = 0
	do i = 1,nj
	   do k = 1,r
	      M(k,xxsub(i)) = M(k,xxsub(i))+V(k,i)*c(i)
	   enddo
           
	enddo
        !part 1 end
ccccccccccccccccccccccccccccccccccccccccccccccccccc
        enddo
        call date_and_time(date,time,zone,values2)
        time1=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_our1         = ',time1/num

        call date_and_time(date,time,zone,values1)
        do j = 1,num

ccccccccccccccccccccccccccccccccccccccccccccccccccc
        !part 2
        
	do i = 1,r
	   in1 = reshape(M(i,:),(/ns,ns,ns/))
	   call dfftw_execute_dft(plan, in1, out1)
	   N(:,:,:,i) = out1
	enddo
        !part 2 end
ccccccccccccccccccccccccccccccccccccccccccccccccccc

        enddo
        call date_and_time(date,time,zone,values2)
        time1=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_our2         = ',time1/num
       
        mm=floor(ns/2.0+0.6)
        !allocate(NN(mm,mm,r))
        !allocate(NN(mm,ns,r))
        !allocate(NN1(ns,mm,r))
        !call date_and_time(date,time,zone,values1)
        !do j = 1,num
        !NN=N(1:mm,1:mm,r)
        !N(1:mm,1:mm,r)=N(mm+1:ns,mm+1:ns,r)
        !N(mm+1:ns,mm+1:ns,r)=NN
        !NN=N(1:mm,mm+1:ns,r)
        !N(1:mm,mm+1:ns,r)=N(mm+1:ns,1:mm,r)
        !N(mm+1:ns,1:mm,r)=NN
        
        !enddo
        !call date_and_time(date,time,zone,values2)
        !time1=sum((values2(5:8)-values1(5:8))*arr)
        !print *,' T_our3         = ',time1/num
       
        call date_and_time(date,time,zone,values1)
        do j = 1,num

cccccccccccccccccccccccccccccccccccccccccccccccccc
          !part 3
           !NN2(1:mm,1:mm,:)=N(1:mm,1:mm,:)
           !NN2(1:mm,mm+1:ns,:)=N(1:mm,ns+mm+1:2*ns,:)
           !NN2(mm+1:ns,1:mm,:)=N(ns+mm+1:2*ns,1:mm,:)
           !NN2(mm+1:ns,mm+1:ns,:)=N(ns+mm+1:2*ns,ns+mm+1:2*ns,:)
           NN1=reshape(N,(/ns*ns*ns,r/))
           !do i = 1,ns*ns
              !S(i)=dot_product(U(i,:),NN1(i,:))
           !enddo
	   S = sum(U*NN1,2)
           !do fftshift
           SS=reshape(S,(/ns,ns,ns/))
           SS1=SS(1:mm,1:mm,1:mm)
           SS(1:mm,1:mm,1:mm)=SS(mm+1:ns,mm+1:ns,mm+1:ns)
           SS(mm+1:ns,mm+1:ns,mm+1:ns)=SS1
           SS1=SS(1:mm,mm+1:ns,1:mm)
           SS(1:mm,mm+1:ns,1:mm)=SS(mm+1:ns,1:mm,mm+1:ns)
           SS(mm+1:ns,1:mm,mm+1:ns)=SS1
           SS1=SS(mm+1:ns,1:mm,1:mm)
           SS(mm+1:ns,1:mm,1:mm)=SS(1:mm,mm+1:ns,mm+1:ns)
           SS(1:mm,mm+1:ns,mm+1:ns)=SS1
           SS1=SS(mm+1:ns,mm+1:ns,1:mm)
           SS(mm+1:ns,mm+1:ns,1:mm)=SS(1:mm,1:mm,mm+1:ns)
           SS(1:mm,1:mm,mm+1:ns)=SS1
           S=reshape(SS,(/ns*ns*ns/))
           !part 3 end
ccccccccccccccccccccccccccccccccccccccccccccccccc          


        enddo
        call date_and_time(date,time,zone,values2)
        time1=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_our3         = ',time1/num
	end subroutine
