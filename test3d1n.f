	program test
	implicit none
        
        integer r,mx
        parameter (mx=500)
        integer ms
        parameter (ms=12)
        integer nj
        parameter (nj=ms*ms*ms)
        integer i,iflag,xsub(nj,3),ier,num,j,xxsub(nj)
        integer n1,n2,mt,k1,k2,mm,k3,rt
        real*16 begin1,end1
        integer*8  time_begin,time_end,countrage,countmax

        real*16 time1,time2
        real*16 arr(4)
        real*8 pi,xj(nj),yj(nj),zj(nj)
        parameter (pi=3.141592653589793238462643383279502884197d0)
        complex*16 U(ms*ms*ms,mx),V(mx,nj),cj(nj),S(ms*ms*ms),re(nj)
        complex*16 fk(ms,ms,ms),fk1(ms*ms*ms),uu(ms,ms,ms,mx)
        complex*16,allocatable :: NN(:,:,:,:)
        real*8 x(nj,3),eps,error
        double complex in1, out1
        dimension in1(ms,ms,ms), out1(ms,ms,ms)
	integer*16 :: plan
        integer FFTW_FORWARD,FFTW_MEASURE
        parameter (FFTW_FORWARD=-1)
        parameter (FFTW_MEASURE=0)
    
        character*8 date
        character*10 time
        character*5 zone 
        integer*4 values1(8),values2(8)
        
        arr(1)=3600
        arr(2)=60
        arr(3)=1
        arr(4)=0.001


        iflag=-1
        eps=1E-12
        num=100
        rt=500





        print *,'start 3D type 1 testing:','nj  =',nj,'ms  =',ms
        print *,'eps             =',eps

        do k3 = -ms/2,(ms-1)/2
	 do k1 = -ms/2, (ms-1)/2
	   do k2 = -ms/2, (ms-1)/2
	      j =  (k1+ms/2+1) + (k2+ms/2)*ms + (k3+ms/2)*ms*ms
	      xj(j) = pi*dcos(-pi*k1/ms)
	      yj(j) = pi*dcos(-pi*k2/ms)
	      zj(j) = pi*dcos(-pi*k3/ms)
	      cj(j) = dcmplx(dsin(pi*j/ms),dcos(pi*j/ms))
	   enddo
	 enddo
        enddo
        x(:,1)=xj
        x(:,2)=yj
        x(:,3)=zj
        call nufft3dI(nj,x,iflag,ms,rt,eps,U,V,xxsub,r)
        print *,'r=',r


        end program
