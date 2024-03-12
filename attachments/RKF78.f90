
!应用主程序示例
program RKF78_test
implicit none
real*8 t,h,T0,Ts,x(6)
integer tsnum
real*8 tpi

tpi=8d0*datan(1d0)

Ts=100d0*tpi
t=0d0
h=Ts/1d3

do while(t<Ts)
	if ( t+h>Ts ) h=Ts-t
	call RKF78(h,t,x,6,1d-14)
enddo

stop
end



!右函数，变量维数需与RKF78中ne相同

	subroutine yhc(t,x,y)
	implicit real*8(a-h,k-z)
	dimension x(6),y(6)

	r=dsqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))

	y(1)=x(4)
	y(2)=x(5)
	y(3)=x(6)
	y(4)=-x(1)/r**3
	y(5)=-x(2)/r**3
	y(6)=-x(3)/r**3

	return
	end




!RKF78积分器

	SUBROUTINE RKF78(H,T,X,ne,ERR)
	IMPLICIT real*8(a-h,k-z)
	integer ne	!Number of equations
	DIMENSION X(ne),X0(ne),X1(ne),Y0(ne),Y1(ne),Y2(ne),Y3(ne),Y4(ne),Y5(ne),Y6(ne),Y7(ne),Y8(ne),Y9(ne),Y10(ne),Y11(ne),Y12(ne)
	T0=T
	DO 20 I=1,ne
20	X0(I)=X(I)
17	DO 16 I=1,ne
	X(I)=X0(I)
16	X1(I)=X(I)
	T=T0
	T1=T
	CALL YHC(T1,X1,Y0)
	DO 1 I=1,ne
1	X1(I)=X(I)+H*2D0/27D0*Y0(I)
	T1=T+H*2D0/27D0
	CALL YHC(T1,X1,Y1)
	DO 2 I=1,ne
2	X1(I)=X(I)+H*(Y0(I)+3D0*Y1(I))/36D0
	T1=T+H*1D0/9D0
	CALL YHC(T1,X1,Y2)
	DO 3 I=1,ne
3	X1(I)=X(I)+H*(Y0(I)+3D0*Y2(I))/24D0
	T1=T+H*1D0/6D0
	CALL YHC(T1,X1,Y3)
	DO 4 I=1,ne
4	X1(I)=X(I)+H*(Y0(I)*20D0+(-Y2(I)+Y3(I))*75D0)/48D0
	T1=T+H*5D0/12D0
	CALL YHC(T1,X1,Y4)
	DO 5 I=1,ne
5	X1(I)=X(I)+H*(Y0(I)+Y3(I)*5D0+Y4(I)*4D0)/20D0
	T1=T+H*1D0/2D0
	CALL YHC(T1,X1,Y5)
	DO 6 I=1,ne
6	X1(I)=X(I)+H*(-Y0(I)*25D0+Y3(I)*125D0-Y4(I)*260D0+Y5(I)*250D0)/108D0
	T1=T+H*5D0/6D0
	CALL YHC(T1,X1,Y6)
	DO 7 I=1,ne
7	X1(I)=X(I)+H*(Y0(I)*93D0+Y4(I)*244D0-Y5(I)*200D0+Y6(I)*13D0)/900D0
	T1=T+H*1D0/6D0
	CALL YHC(T1,X1,Y7)
	DO 8 I=1,ne
8	X1(I)=X(I)+H*(Y0(I)*180D0-Y3(I)*795D0+Y4(I)*1408D0-Y5(I)*1070D0+Y6(I)*67D0+Y7(I)*270D0)/90D0
	T1=T+H*2D0/3D0
	CALL YHC(T1,X1,Y8)
	DO 9 I=1,ne
9	X1(I)=X(I)+H*(-Y0(I)*455D0+Y3(I)*115D0-Y4(I)*3904D0+Y5(I)*3110D0-Y6(I)*171D0+Y7(I)*1530D0-Y8(I)*45D0)/540D0
	T1=T+H*1D0/3D0
	CALL YHC(T1,X1,Y9)
	DO 10 I=1,ne
10	X1(I)=X(I)+H*(Y0(I)*2383D0-Y3(I)*8525D0+Y4(I)*17984D0-Y5(I)*15050D0+Y6(I)*2133D0+Y7(I)*2250D0+Y8(I)*1125D0+Y9(I)*1800D0)/4100D0
	T1=T+H
	CALL YHC(T1,X1,Y10)
	DO 11 I=1,ne
11	X1(I)=X(I)+H*(Y0(I)*60D0-Y5(I)*600D0-Y6(I)*60D0+(Y8(I)-Y7(I)+2D0*Y9(I))*300D0)/4100D0
	T1=T
	CALL YHC(T1,X1,Y11)
	DO 12 I=1,ne
12	X1(I)=X(I)+H*(-Y0(I)*1777D0-Y3(I)*8525D0+Y4(I)*17984D0-Y5(I)*14450D0+Y6(I)*2193D0+Y7(I)*2550D0+Y8(I)*825D0+Y9(I)*1200D0+Y11(I)*4100D0)/4100D0
	T1=T+H
	CALL YHC(T1,X1,Y12)
	D=0D0
	DO 14 I=1,ne
	X(I)=X(I)+H*(Y5(I)*272D0+(Y6(I)+Y7(I))*216D0+(Y8(I)+Y9(I))*27D0+(Y11(I)+Y12(I))*41D0)/840D0
14	D=D+DABS((Y0(I)+Y10(I)-Y11(I)-Y12(I))*H*41D0/840D0)
	T=T+H
	IF(D.GT.ERR) THEN
	H=H/2D0
	GOTO 17
	END IF
	IF(D.LE.ERR*1D-4.and.h.le.5d-3) H=H*2D0
	RETURN
	END


