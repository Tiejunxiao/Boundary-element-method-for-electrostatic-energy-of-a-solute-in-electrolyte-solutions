CC*************************************************************************************************************************
CC*   Copyright (c) 2026 [Xueyu Song and Tiejun Xiao]. 
CC*  This code is distributed under the terms of the Creative Commons Attribution 4.0 International  License (CC BY 4.0).
CC*  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
CC*
CC*  Author: Xueyu Song
CC*  Email:   [xsong@iastate.edu]
CC*  Author: Tiejun Xiao
CC*  Email:   [tjxiao@gznc.edu.cn]
CC*  Date:    2026-10-04
CC*************************************************************************************************************************
CC*  This code uses the boundary element method(BEM) to solve the  Debye-Hückel(DH)  equation
CC*    and evaluate the electrostatic energy of a tetraatomic solute in an electrolyte solution.
CC*  ref :[1] A. H. Juffer, E. F. F. Botta, B. A. M. Van Keulen, A. Van der Ploeg, and H. J. C. Berendsen, The
CC*     electric potential of a macromolecule in a solvent: a fundamental approach, J. Comput. Phys. 97, 144 (1991)
CC*   [2] K. E. Atkinson, The numerical solution of integral equations of the second kind (Cambridge University Press, 1997)
CC*   
CC*  The DH parameter of the electrolyte solution is kDH=0.805/Angstrom.
CC*  This code computes  the  electrostatic energy  ue for the DH mode with inverse Debye length kDH. 
CC*************************************************************************************************************************
      IMPLICIT none
      integer node,nedg,nele
      integer mxface,mxvert,machin,nc,M
      parameter(node=618,nedg=1848,nele=1232)
      PARAMETER (MXFACE =nele, MXVERT = 2*(MXFACE+1), MACHIN=MXVERT)
      parameter(nc=4,M=1)
CC*  nc is the numberof atoms of the solute 
      character option*6
      integer IFACE(7,MXFACE),IPIVTN(MXVERT*2)
      integer LEVEL(MXVERT),ipair(2,nedg)
      complex VRINT(6)
      real  RADIUS(MXFACE),CENTER(3,MXFACE), VERTEX(3,MXVERT)
      real  vernor(3,mxvert),Rj(3),W(3),c343(4,nc),P(3)
      real xi(3),xj(3),ai(3),aj(3),xh(3),ah(3),pn(3),h_min,h_max,ue
      complex  KERMAT(MXVERT*2,MXVERT*2), RHS(MXVERT*2,1)
      complex xk,eps,solv,U
      real sum,dd,ds,zero,one,two,four,pi,twopi,fourpi,eps1,dc,E_trans
      integer i,mi,mj,mk,iadd,k,n1,n2,n3,l,mxlevl,ninteg,m1,m2,nls
      integer ngt,iflag,ii,nface,nvert,iw,j,jz,ind,jx,ix,info
      common/block2/Rj
      common/norm/pn,xk,eps

C     **** INITIALIZATION ****
C
CC*  En_trans= 1389.35 is an energy  coefficient for electrostatic interactions.
CC*  i.e., the electrostatic energy of two elementary charges separated by 1 Angstrom is about 1389.35  kJ/mol.
CC*  1 e_0^2/(4 *PI * Epsilon_0* 1 A) =1389.35 kJ/mol
CC*  e_0 is the elementary charge
CC*  Epsilon_0 is the permittivity of the vacuum.
CC*  eps1= 78.0 is the relative dielectric constant of the aqueous electrolyte solution. 

      E_trans= 1389.35 
      eps1=78.0   
      zero=0.0
      one=1.0
      two=2.0
      four=4.0
      PI = 3.141592653
      TWOPI = TWO*PI
      FOURPI = FOUR*PI
C
C     **** INPUT OF PROBLEM PARAMETERS ****

       open(unit=30,file='FA.crg',status='old')
       open(unit=40,file='FA.vert',status='old')
       open(unit=50,file='FA.face',status='old')

      open(unit=20,file='ue-DH.dat',status='unknown')

        do 30 i=1,nc
          read(30,*)c343(1,i),c343(2,i),c343(3,i),c343(4,i)
30      continue

cccccccccccc   define the triangulation of the molecular surface cccc 

      read(40,*)
      read(40,*)
      read(40,*)
      do i=1,node
        read(40,*)vertex(1,i),vertex(2,i),vertex(3,i),
     *    vernor(1,i),vernor(2,i),vernor(3,i),mi,mj,mk
      enddo

      write(*,*)'Reading triangulation files0...'
      iadd=0
      read(50,*)
      read(50,*)
      read(50,*)
C   read in three indices of the corner vertices of face k
      do 592 k=1,nele
        read(50,*)n1,n2,n3,mi,mj
C   store them in iface
        iface(1,k)=n1
        iface(2,k)=n2
        iface(3,k)=n3

C   VERTEX #4,#5,#6
       do 542 L=4,6
        if(L.eq.4) then
          m1=n1
          m2=n2
        else if(L.eq.5) then
          m1=n2
          m2=n3
        else if(L.eq.6) then
          m1=n1
          m2=n3
        end if

        if(m1.gt.m2) then
          nls=m2
          ngt=m1
        else
          nls=m1
          ngt=m2
        end if

        iflag=0
        do ii=1,iadd
          if(nls.eq.ipair(1,ii)
     *        .and. ngt.eq.ipair(2,ii)) then
             iflag=1
             iface(L,k)=ii+node
          endif
        end do
        if(iflag.eq.0) then
          iadd=iadd+1
          ind=iadd+node
          iface(L,k)=ind
          xi(1)=vertex(1,m1)
          xi(2)=vertex(2,m1)
          xi(3)=vertex(3,m1)
          xj(1)=vertex(1,m2)
          xj(2)=vertex(2,m2)
          xj(3)=vertex(3,m2)

          ai(1)=vernor(1,m1)
          ai(2)=vernor(2,m1)
          ai(3)=vernor(3,m1)
          aj(1)=vernor(1,m2)
          aj(2)=vernor(2,m2)
          aj(3)=vernor(3,m2)

          ipair(1,iadd)=nls
          ipair(2,iadd)=ngt

          call param(xi,xj,ai,aj,xh,ah)

          vertex(1,ind)=xh(1)
          vertex(2,ind)=xh(2)
          vertex(3,ind)=xh(3)

          vernor(1,ind)=ah(1)
          vernor(2,ind)=ah(2)
          vernor(3,ind)=ah(3)

        end if

542    continue
592   continue
        nface=mxface
        nvert=mxvert

         CALL CIRCLE(IFACE,NFACE,VERTEX,RADIUS,CENTER)
         DO 15 I=1,NFACE
15        RADIUS(I) = SQRT(RADIUS(I))
          CALL DIAMTR(IFACE,NFACE,VERTEX,H_MIN,H_MAX)

cccccccccccccc end of triangulation ccccccccccccccccccccccccccccc

         option='CURVED'
        mxlevl=2
        ninteg=10

C***  eps=1.0  represents  the ratio of dielectric constant of the solvent and the solute
C***  xk=0.805/Angstrom  is the inverse  Debye length of the electrolyte solution           
          eps=1.0
          xk=cmplx(0.805,0.00)

      do j=1,nvert
        sum=0.0
        do  jx=1,nc
               dd=sqrt((c343(1,jx)-vertex(1,j))**2+
     *                 (c343(2,jx)-vertex(2,j))**2+
     *                 (c343(3,jx)-vertex(3,j))**2)

            sum=sum+c343(4,jx)/dd
          enddo
        RHS(j,1)=sum/fourpi/eps1
          enddo

      do 220 j=1,nvert
        sum=zero
        do 221 jx=1,nc
               dd=sqrt((c343(1,jx)-vertex(1,j))**2+
     *                 (c343(2,jx)-vertex(2,j))**2+
     *                 (c343(3,jx)-vertex(3,j))**2)

            dc=c343(1,jx)*vernor(1,j)+
     *         c343(2,jx)*vernor(2,j)+
     *         c343(3,jx)*vernor(3,j)
            ds=vertex(1,j)*vernor(1,j)+
     *         vertex(2,j)*vernor(2,j)+
     *         vertex(3,j)*vernor(3,j)
            sum=sum+c343(4,jx)*(dc-ds)/dd**3
221     continue
        RHS(j+nvert,1)=sum/fourpi/eps1
220   continue
             
      DO 20 I=1,NVERT*2
      DO 20 J=1,NVERT*2
20     KERMAT(I,J) = cmplx(0.0,0.0)

C       kernel 00 part
      DO 80 K=1,NFACE
        CALL SETLEV(K,IFACE,VERTEX,NVERT,CENTER,RADIUS,
     *       LEVEL,MXLEVL,H_MAX)
        DO 70 I=1,NVERT
          DO 60 L=1,3
60          P(L) = VERTEX(L,I)
          IF(LEVEL(I) .GE. 0) THEN
C             P IS OUTSIDE FACE #K.
              CALL INTKL2m(K,P,VRINT,OPTION,LEVEL(I),
     *               IFACE,VERTEX,'DOUBLE',xk,eps)
          ELSE
C             P IS INSIDE FACE #K.
              CALL INTKL4m(K,I,VRINT,OPTION,NINTEG,'DOUBLE',
     *                        IFACE,VERTEX,xk,eps)
          END IF
          DO  J=1,6
            L = IFACE(J,K)
            KERMAT(I,L) = KERMAT(I,L)-VRINT(J)/fourpi
           enddo
70        continue

C       kernel 01 part
        DO 81 I=1,NVERT
          DO 61 L=1,3
61          P(L) = VERTEX(L,I)
          IF(LEVEL(I) .GE. 0) THEN
C             P IS OUTSIDE FACE #K.
              CALL INTKL2m(K,P,VRINT,OPTION,LEVEL(I),
     *           IFACE,VERTEX,'SINGLE',xk,eps)
          ELSE
C             P IS INSIDE FACE #K.
              CALL INTKL4m(K,I,VRINT,OPTION,NINTEG,'SINGLE',
     *                        IFACE,VERTEX,xk,eps)
          END IF
          DO  J=1,6
            L = IFACE(J,K)+nvert
            KERMAT(I,L) = KERMAT(I,L)-VRINT(J)/fourpi
               enddo
81        continue

c            10 part
        DO 72 I=1,NVERT
          DO 62 L=1,3
            P(L) = VERTEX(L,I)
62          pn(L) =vernor(L,i)
          IF(LEVEL(I) .GE. 0) THEN
C             P IS OUTSIDE FACE #K.
              CALL INTKLm(K,P,pn,VRINT,OPTION,LEVEL(I),
     *               IFACE,VERTEX,'DOUBLE',xk,eps)
          ELSE
C             P IS INSIDE FACE #K.
              CALL INTKLmm(K,I,pn,VRINT,OPTION,NINTEG,'DOUBLE',
     *               IFACE,VERTEX,xk,eps)
          END IF
          DO  J=1,6
            L = IFACE(J,K)
         KERMAT(I+nvert,L) = KERMAT(I+nvert,L)-VRINT(j)/fourpi
          enddo
72       continue
                          
cccccccccc         11 part
        DO 73 I=1,NVERT
          DO 63 L=1,3
            P(L) = VERTEX(L,I)
63          pn(L) =vernor(L,i)
          IF(LEVEL(I) .GE. 0) THEN
C             P IS OUTSIDE FACE #K.
              CALL INTKL1m(K,P,VRINT,OPTION,LEVEL(I),
     *                        IFACE,VERTEX)
          ELSE
C             P IS INSIDE FACE #K.
              CALL INTKL3m(K,I,VRINT,OPTION,NINTEG,
     *                        IFACE,VERTEX)
          END IF
          DO  J=1,6
            L = IFACE(J,K)+nvert
            KERMAT(I+nvert,L) = KERMAT(I+nvert,L)-VRINT(J)/fourpi
            enddo
73        continue
80          continue

        DO 100 I=1,NVERT
100     KERMAT(I,I) = KERMAT(I,I) +0.5*(one+eps)

        DO 101 I=nvert+1,NVERT*2
101     KERMAT(I,I) = KERMAT(I,I) +0.5*(one+one/eps)

         
c       write(*,*)'first stage is ok'
C     SOLVE LINEAR SYSTEM.
         call cgetrf(nvert*2,nvert*2,kermat,nvert*2,ipivtn,info)
           write(*,*)info
         if(info.ne.zero) then
            write(20,*)'the matrix is singular'
            stop 
           endif
                write(*,*)'pass5'
      call cgetrs('N',nvert*2,1,kermat,nvert*2,
     *            ipivtn,RHS,nvert*2,info)

         if(info.ne.zero) then
            write(20,*)'the solution is singular'
            stop 
           endif

c        calculate the f(a,a)
CC     solv saves the electrostatic energy from DH theory

         solv=cmplx(0.0,0.0)
        do 900 ix=1,nc
             P(1)=c343(1,ix)
             P(2)=c343(2,ix)
             P(3)=c343(3,ix)
       CALL LEVDIR(P,IFACE,NFACE,CENTER,LEVEL,MXLEVL,H_MAX)
      U = cmplx(0.0,0.0)
         DO 800 K=1,NFACE
           CALL INTKL2m(K,P,VRINT,OPTION,LEVEL(k),
     *            IFACE,VERTEX,'DOUBLE',xk,eps)
            DO 840 J=1,6
               L = IFACE(J,K)
840            U = U + VRINT(J)*RHS(L,1)
800      CONTINUE
         DO 801 K=1,NFACE
           CALL INTKL2m(K,P,VRINT,OPTION,LEVEL(k),
     *            IFACE,VERTEX,'SINGLE',xk,eps)
            DO 841 J=1,6
               L = IFACE(J,K)+nvert
841           U = U + VRINT(J)*RHS(L,1)
801     CONTINUE
           solv=solv+U*c343(4,ix)
900   continue

           ue= Real(solv/2.*E_trans)
         write(20,1002) ue
         write(*,1002) ue
              
1002  format('The electrostatic energy from the DH theory is ', f12.8,
     * ' kJ/mol.') 

      END

      function kernel(p,q)
C     ------------------------------------------------
      IMPLICIT none
      complex KERNEL
      real P(3),Q(3),pn(3),qp(3),one,dd,up
      complex xk,eps
      common/norm/pn,xk,eps
      one=1.0
      dd=(p(1)-q(1))**2+(p(2)-q(2))**2+(p(3)-q(3))**2
      dd=sqrt(dd)
      qp(1)=q(1)-p(1)
      qp(2)=q(2)-p(2)
      qp(3)=q(3)-p(3)
      up=pn(1)*qp(1)+pn(2)*qp(2)+pn(3)*qp(3)
      kernel =up/dd**3*(one-(one+xk*dd)*exp(-xk*dd)/eps)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LEVDIR(P,IFACE,NFACE,CENTER,LEVEL,MAXLEV,H)
      IMPLICIT none 
      
      integer IFACE(7,*), LEVEL(*),K,L,NFACE,MAXLEV
      real CENTER(3,*),P(3),R,H
C
      DO 100 K=1,NFACE
	R = SQRT((P(1) - CENTER(1,K))**2
     *	       + (P(2) - CENTER(2,K))**2
     *	       + (P(3) - CENTER(3,K))**2)
	L = INT(R/H)
	LEVEL(K) = MAX( 0, MAXLEV - L )
100	CONTINUE
      RETURN
      END

       subroutine param(xi,xj,ai,aj,xh,ah)
      implicit none
      real xi(3),xj(3),ai(3),aj(3),xh(3),ah(3)
      real dx1(3),dx2(3),am(3)
      real zero,one,two,three,four,five,six,eight
      parameter(zero=0.0,one=1.0,two=2.0,three=3.0)
      parameter(four=4.0,five=5.0,six=6.0)
      parameter(eight=four*two)
      real aij,xjiai,xjiaj,alpha,beta,amag,d1,d2,d12
      integer l

      aij=ai(1)*aj(1)+ai(2)*aj(2)+ai(3)*aj(3)
      xjiai=xj(1)*ai(1)+xj(2)*ai(2)+xj(3)*ai(3)
     *          -xi(1)*ai(1)-xi(2)*ai(2)-xi(3)*ai(3)
      xjiaj=xj(1)*aj(1)+xj(2)*aj(2)+xj(3)*aj(3)
     *          -xi(1)*aj(1)-xi(2)*aj(2)-xi(3)*aj(3)
      alpha=(six*xjiai+three*xjiaj*aij)/(four-aij)
      beta=-(three*xjiaj+alpha*aij)/two

      do L=1,3
        dx1(L)=xj(L)-xi(L)-
     *          one/(four*three)*(beta*aj(L)-alpha*ai(L))
        dx2(L)=alpha*ai(L)+beta*aj(L)
        xh(L)=(xi(L)+xj(L))/two-dx2(L)/eight
      end do

      d1=sqrt(dx1(1)**2+dx1(2)**2+dx1(3)**2)
      d2=sqrt(dx2(1)**2+dx2(2)**2+dx2(3)**2)
      d12=dx1(1)*dx2(1)+dx1(2)*dx2(2)+dx1(3)*dx2(3)
      amag=zero
      do L=1,3
        ah(L)=(dx2(L)-(d12/d1**2)*dx1(L))/d1
        amag=amag+ah(L)**2
      end do
      amag=sqrt(amag)

      do L=1,3
        am(L)=ai(L)+aj(L)
        ah(L)=ah(L)/amag
      end do
      if(am(1)*ah(1)+am(2)*ah(2)+am(3)*ah(3).lt.zero) then
        do L=1,3
          ah(L)=-ah(L)
        end do
      end if

      return
      end

      SUBROUTINE MIDPT(IDX,V,IVX,I1,I2,LED,IFACE,VERTEX,INDXV)
C     ----------------
C
C     FIND THE MIDPOINT 'V' BETWEEN TWO GIVEN NODES, OF INDICES
C     #I1 AND #I2 WITHIN FACE #IDX OF 'IFACE'.  ALSO SET 'IVX' TO
C     INDICATE WHETHER 'V' IS AT A VERTEX, EDGE, OR FACE OF THE
C     ORIGINAL SURFACE S.  LED=0 MEANS THAT 'V' IS NOT ON AN EDGE
C     OF S NOR IS IT A VERTEX OF S; AND THUS IVX=0 IN THIS CASE.
C
      IMPLICIT real (A-H,O-Z)
      DIMENSION V(3), IFACE(7,*), VERTEX(3,*), INDXV(*)
      DATA TWO/2.0/

      return
      end

      SUBROUTINE CIRCLE(IFACE,NFACE,VERTEX,RADSQR,CENTER)
      IMPLICIT real (A-H,O-Z)
      PARAMETER (ZERO=0.0, ONE=1.0, TWO=2.0)
      DIMENSION IFACE(7,*), VERTEX(3,*), RADSQR(*), CENTER(3,*),
     *          A(3,3), B(3), DIFF(3), IV(3), SUM(3), IPIVOT(3),
     *          WORKING(3)
C
C     CONSIDER FACE #I.
      DO 100 I=1,NFACE
C       INITIALIZE.
        DO 10 J=1,3
          A(1,J) = ONE
10        IV(J) = IFACE(J,I)
        B(1) = ONE
        K = IV(1)
C       INITIALIZE FOR COMPUTATION OF EQUATIONS 2 AND 3.
        DO 50 L=2,3
          M = IV(L)
          DO 20 J=1,3
            DIFF(J) = VERTEX(J,M) - VERTEX(J,K)
20          SUM(J) = (VERTEX(J,M) + VERTEX(J,K))/TWO
          B(L) = DIFF(1)*SUM(1) + DIFF(2)*SUM(2) + DIFF(3)*SUM(3)
C         SET UP COEFFICIENTS OF ROW #L OF A.
          DO 30 J=1,3
            N = IV(J)
30          A(L,J) = DIFF(1)*VERTEX(1,N) + DIFF(2)*VERTEX(2,N)
     *             + DIFF(3)*VERTEX(3,N)
50        CONTINUE
C       SOLVE FOR BARYMETRIC COORDINATES.
c        CALL LINSYS(A,A,3,B,B,1,RCOND,IPIVOT,3,WORKING,ERRMAX)
        call SGETRF( 3,3, A, 3, IPIVOT, INFO )
        call SGETRS( 'N', 3, 1, A, 3, IPIVOT, B, 3, INFO )
C       OBTAIN RECTANGULAR COORDINATES OF CENTER.
        DO 60 J=1,3
60        CENTER(J,I) = B(1)*VERTEX(J,IV(1)) + B(2)*VERTEX(J,IV(2))
     *                + B(3)*VERTEX(J,IV(3))
C       OBTAIN SQUARE OF RADIUS.
        SM = ZERO
        K = IV(1)
        DO 70 J=1,3
70        SM = SM + (CENTER(J,I) - VERTEX(J,K))**2
        RADSQR(I) = SM
100     CONTINUE
      RETURN
      END

      SUBROUTINE DIAMTR(IFACE,NFACE,VERTEX,H_MIN,H_MAX)
C     SUBPROGRAMS CALLED: NONE
      IMPLICIT real (A-H,O-Z)
      PARAMETER (ZERO=0.0, THOUSN=1000.0)
      DIMENSION IFACE(7,*), VERTEX(3,*), D(3), JP1(3)
      DATA JP1/2,3,1/
      H_MAX = ZERO
      H_MIN = THOUSN
      DO 20 I=1,NFACE
        DO 10 J=1,3
          J1 = IFACE(J,I)
          J2 = IFACE(JP1(J),I)
          D(J) = SQRT((VERTEX(1,J2) - VERTEX(1,J1))**2
     *              + (VERTEX(2,J2) - VERTEX(2,J1))**2
     *              + (VERTEX(3,J2) - VERTEX(3,J1))**2)
10        CONTINUE
        DMAX = MAX(D(1),D(2),D(3))
        H_MIN = MIN(H_MIN,DMAX)
        H_MAX = MAX(H_MAX,DMAX)
20      CONTINUE
      RETURN
      END

      SUBROUTINE SETLEV(K,IFACE,VERTEX,NVERT,CENTER,RADIUS,
     *                  LEVEL,MAXLEV,H)
      IMPLICIT real         (A-H,O-Z)
      DIMENSION IFACE(7,*), VERTEX(3,*), CENTER(3,*), RADIUS(*),
     *          LEVEL(*)
C
      DO 20 I=1,NVERT
20      LEVEL(I) = 0
      DO 30 I=1,6
30      LEVEL(IFACE(I,K)) = -1
      DO 40 I=1,NVERT
        IF(LEVEL(I) .NE. -1) THEN
            R = SQRT((VERTEX(1,I) - CENTER(1,K))**2
     *             + (VERTEX(2,I) - CENTER(2,K))**2
     *             + (VERTEX(3,I) - CENTER(3,K))**2)
            L = MAX( INT((R - RADIUS(K))/H), 0 )
            LEVEL(I) = MAX( 0, MAXLEV - L )
        END IF
40      CONTINUE
      RETURN
      END

        SUBROUTINE INTKL2m(IDX,P,VRINT,OPTION,LEVEL,
     *                     IFACE,VERTEX,LAYER,xk,eps)
      IMPLICIT none
      CHARACTER LAYER*6, LYR*6, OPTION*6, IOP*6
      complex VRINT,xxk,xeps,xk,eps
      integer IFACE,i,k,j,idx,level
      real VERTEX,VRTX,PDUP,P,VS,VT,VN,AREA
      DIMENSION IFACE(7,*), VERTEX(3,*), VRINT(6), VRTX(3,6), PDUP(3),
     *          P(3), VS(3), VT(3)
      COMMON/PFNKL2m/VRTX,AREA,PDUP,VN(3),IOP,LYR,xxk,xeps
      EXTERNAL FCNKL2m
C
C     INITIALIZE COMMON FOR DEFINITION OF FCNKL2.
      IOP = OPTION
      LYR = LAYER
      xxk=xk
      xeps=eps
      DO 10 I=1,3
10      PDUP(I) = P(I)
      DO 20 J=1,6
        K = IFACE(J,IDX)
        DO 20 I=1,3
20        VRTX(I,J) = VERTEX(I,K)
      IF(IOP .EQ. 'PLANAR') THEN
          DO 30 I=1,3
            VT(I) = VRTX(I,2) - VRTX(I,1)
30          VS(I) = VRTX(I,3) - VRTX(I,1)
          IF(LAYER .EQ. 'SINGLE') THEN
C             CALCULATE AREA FOR PLANAR TRIANGLE.
              AREA = SQRT((VS(1)*VT(2) - VS(2)*VT(1))**2
     *                  + (VS(1)*VT(3) - VS(3)*VT(1))**2
     *                  + (VS(2)*VT(3) - VS(3)*VT(2))**2)
          ELSE
C             COMPUTE AN ORTHOGONAL VECTOR.
              VN(1) = VT(2)*VS(3) - VT(3)*VS(2)
              VN(2) = VT(3)*VS(1) - VT(1)*VS(3)
              VN(3) = VT(1)*VS(2) - VT(2)*VS(1)
          END IF
      END IF
C
C     PERFORM INTEGRATION.
      CALL SMPLX6m(VRINT,FCNKL2m,LEVEL)
      RETURN
      END

      SUBROUTINE FCNKL2m(S,T,F)
C     -----------------
      IMPLICIT none
      CHARACTER LAYER*6, OPTION*6
      real zero,one,two,three,four
      PARAMETER (ZERO=0.0, ONE=1.0, TWO=2.0, THREE=3.0,
     *           FOUR=4.0)
      complex F,FK,xk,eps
      real v,p,q,c,dqs,dqt,dot,dist,vn,area,s,t,spt,CCUTRD,u,sum
      integer i,j,k
      DIMENSION V(3,6), P(3), Q(3), F(6), C(6), DQS(3), DQT(3)
      COMMON/PFNKL2m/V,AREA,P,VN(3),OPTION,LAYER,xk,eps
C
C     THE FOLLOWING CONSTANT IS TO BE A SMALL NUMBER, OF THE ORDER
C     OF THE UNIT ROUND, SAY 1000 TIMES THE UNIT ROUND.
c      CCUTRD = 1000.0*D1MACH(3)
      CCUTRD = 1000.0*1.81899e-32
C
C     CALCULATE CARDINAL FUNCTIONS C(I,S,T).
      SPT = S + T
      U = ONE - SPT
      C(1) = U*(TWO*U - ONE)
      C(2) = T*(TWO*T - ONE)
      C(3) = S*(TWO*S - ONE)
      C(4) = FOUR*T*U
      C(5) = FOUR*S*T
      C(6) = FOUR*S*U
C
C     CALCULATE Q=Q(S,T).
      DO 20 I=1,3
        SUM = ZERO
        DO 10 J=1,6
10        SUM = SUM + C(J)*V(I,J)
20      Q(I) = SUM
C
C     DO SURFACE ELEMENT AREA DIFFERENTIAL OR ORTHOGONAL VECTOR.
      IF(OPTION .EQ. 'CURVED') THEN
          DO 30 I=1,3
            DQS(I) = (FOUR*SPT - THREE)*V(I,1) + (FOUR*S - ONE)*V(I,3)
     *             + FOUR*(-T*V(I,4) + T*V(I,5)
     *             + (ONE - T - TWO*S)*V(I,6))
30          DQT(I) = (FOUR*SPT - THREE)*V(I,1) + (FOUR*T - ONE)*V(I,2)
     *             + FOUR*((ONE - S - TWO*T)*V(I,4)
     *             + S*(V(I,5) - V(I,6)))
          IF(LAYER .EQ. 'SINGLE') THEN
C             CALCULATE AREA DIFFERENTIAL.
              AREA = SQRT((DQS(1)*DQT(2) - DQS(2)*DQT(1))**2
     *                  + (DQS(1)*DQT(3) - DQS(3)*DQT(1))**2
     *                  + (DQS(2)*DQT(3) - DQS(3)*DQT(2))**2)
          ELSE
C             CALCULATE ORTHOGONAL VECTOR.
              VN(1) = DQT(2)*DQS(3) - DQT(3)*DQS(2)
              VN(2) = DQT(3)*DQS(1) - DQT(1)*DQS(3)
              VN(3) = DQT(1)*DQS(2) - DQT(2)*DQS(1)
          END IF
      END IF
C
C     CALCULATE KERNEL(P,Q).
      DIST = SQRT((P(1) - Q(1))**2 + (P(2) - Q(2))**2
     *          + (P(3) - Q(3))**2)
      IF(DIST .LE. CCUTRD) THEN
C         BECAUSE OF THE NEARNESS OF P TO Q, DEFINE KERNEL(P,Q)=ZERO.
          DO 40 I=1,6
40          F(I) = cmplx(0.0,0.0)
          RETURN
      END IF
      IF(LAYER .EQ. 'SINGLE') THEN
C         CALCULATE THE SINGLE LAYER KERNEL.
          FK = AREA/DIST*(one-exp(-xk*dist))
      ELSE
C         CALCULATE THE DOUBLE LAYER KERNEL.
          DOT = VN(1)*(P(1) - Q(1)) + VN(2)*(P(2) - Q(2))
     *        + VN(3)*(P(3) - Q(3))
          FK = DOT/DIST**3*(eps*(one+xk*dist)*exp(-xk*dist)-one)
      END IF
C
C     CALCULATE INTEGRAND.
      DO 50 I=1,6
50      F(I) = FK*C(I)
      RETURN
      END


      SUBROUTINE INTKL4m(IDXF,IDXP,VRINT,OPTION,NINTEG,LAYER,
     *                  IFACE,VERTEX,xk,eps)
C
      IMPLICIT none
      CHARACTER OPTION*6, IOP*6, LAYER*6, LYR*6
      real two
      PARAMETER (TWO=2.0)
      complex xxk,xeps,xk,eps,VRINT,VITEMP
      integer IFACE,i,j,k,IPOINT,IDXP2,IDXF,IDXP,NINTEG
      real VERTEX,VRTX,P,VS,VT,VN,AREA
      DIMENSION IFACE(7,*), VERTEX(3,*), VRINT(6), VRTX(3,6), P(3),
     *          VS(3), VT(3), VITEMP(6), VN(3)
      COMMON/PFNKL4m/VRTX,AREA,P,IOP,LYR,IPOINT,VN,xxk,xeps
      EXTERNAL FCNKL4m
C
C     INITIALIZE COMMON FOR DEFINITION OF FCNKL4.
      IOP = OPTION
      LYR = LAYER
      xxk=xk
      xeps=eps
      DO 10 I=1,3
10      P(I) = VERTEX(I,IDXP)
      DO 20 J=1,6
        K = IFACE(J,IDXF)
        DO 20 I=1,3
20        VRTX(I,J) = VERTEX(I,K)
      IF(OPTION .EQ. 'PLANAR') THEN
          DO 30 I=1,3
            VT(I) = VRTX(I,2) - VRTX(I,1)
30          VS(I) = VRTX(I,3) - VRTX(I,1)
          IF(LAYER .EQ. 'SINGLE') THEN
C             CALCULATE AREA FOR PLANAR TRIANGLE.
              AREA = SQRT((VS(1)*VT(2) - VS(2)*VT(1))**2
     *                  + (VS(1)*VT(3) - VS(3)*VT(1))**2
     *                  + (VS(2)*VT(3) - VS(3)*VT(2))**2)
          ELSE
C             COMPUTE AN ORTHOGONAL VECTOR.
              VN(1) = VT(2)*VS(3) - VT(3)*VS(2)
              VN(2) = VT(3)*VS(1) - VT(1)*VS(3)
              VN(3) = VT(1)*VS(2) - VT(2)*VS(1)
          END IF
      END IF
C
C     LOCATE INDEX OF P RELATIVE TO THE FACE #IDXF.
      DO 40 I=1,6
        IF(IDXP .EQ. IFACE(I,IDXF)) THEN
            IDXP2 = I
            GO TO 45
        END IF
40      CONTINUE
      PRINT *, ' ERROR IN SUBROUTINE INTKL4.'
      PRINT *, ' THE POINT P IS NOT IN THE GIVEN TRIANGULAR FACE.'
      STOP
C
C     PERFORM INTEGRATION.
45    IF(IDXP2 .LT. 4) THEN
          IPOINT = IDXP2
          CALL INTEGRm(VRINT,FCNKL4m,NINTEG)
          RETURN
      ELSE
          IPOINT = 2*IDXP2 - 4
          CALL INTEGRm(VITEMP,FCNKL4m,NINTEG)
          IPOINT = IPOINT + 1
          CALL INTEGRm(VRINT,FCNKL4m,NINTEG)
          DO 50 J=1,6
50          VRINT(J) = (VRINT(J) + VITEMP(J))/TWO
          RETURN
      END IF
      END

      SUBROUTINE FCNKL4m(X,Y,F)
C
      IMPLICIT none
      CHARACTER OPTION*6, LAYER*6
      real zero,one,two,three,four
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0,
     *           FOUR=4.0D0)
      complex xk,eps,F,FK
      real v,p,q,c,dqs,dqt,vn,area,dot,u,s,t,spt,CCUTRD,x,y,sum,dist
      integer i,j,k,IPOINT
      DIMENSION V(3,6), P(3), Q(3), F(6), C(6), DQS(3), DQT(3), VN(3)
      COMMON/PFNKL4m/V,AREA,P,OPTION,LAYER,IPOINT,VN,xk,eps
C
C     THE FOLLOWING CONSTANT IS TO BE A SMALL NUMBER, OF THE ORDER
C     OF THE UNIT ROUND , SAY 1000 TIMES THE UNIT ROUND.
c      CCUTRD = 1000*D1MACH(3)
        CCUTRD = 1000.0*1.81899e-32
C
C     CALCULATE (S,T) FROM (X,Y).
      GO TO (1,2,3,4,5,6,7,8,9), IPOINT
1     S = Y*X
      T = (ONE - Y)*X
      GO TO 10
2     S = Y*X
      T = ONE - X
      GO TO 10
3     S = ONE - X
      T = (ONE - Y)*X
      GO TO 10
4     S = Y*X
      T = (ONE - X)/TWO + (ONE - Y)*X
      GO TO 10
5     S = (ONE - Y)*X
      T = (ONE - X)/TWO
      GO TO 10
6     S = (ONE - X)/TWO
      T = S + (ONE - Y)*X
      GO TO 10
7     T = (ONE - X)/TWO
      S = T + Y*X
      GO TO 10
8     S = (ONE - X)/TWO
      T = Y*X
      GO TO 10
9     S = (ONE - X)/TWO + Y*X
      T = (ONE - Y)*X
C
C     CALCULATE CARDINAL FUNCTIONS C(*,S,T).
10    SPT = S + T
      U = ONE - SPT
      C(1) = U*(TWO*U - ONE)
      C(2) = T*(TWO*T - ONE)
      C(3) = S*(TWO*S - ONE)
      C(4) = FOUR*T*U
      C(5) = FOUR*S*T
      C(6) = FOUR*S*U
C
C     CALCULATE Q=Q(S,T).
      DO 20 I=1,3
        SUM = ZERO
        DO 15 J=1,6
15        SUM = SUM + C(J)*V(I,J)
20      Q(I) = SUM
C
C     DO SURFACE ELEMENT AREA DIFFERENTIAL OR ORTHOGONAL VECTOR.
      IF(OPTION .EQ. 'CURVED') THEN
          DO 30 I=1,3
            DQS(I) = (FOUR*SPT - THREE)*V(I,1) + (FOUR*S - ONE)*V(I,3)
     *             + FOUR*(-T*V(I,4) + T*V(I,5)
     *             + (ONE - T - TWO*S)*V(I,6))
30          DQT(I) = (FOUR*SPT - THREE)*V(I,1) + (FOUR*T - ONE)*V(I,2)
     *             + FOUR*((ONE - S - TWO*T)*V(I,4)
     *             + S*V(I,5) - S*V(I,6))
          IF(LAYER .EQ. 'SINGLE') THEN
C             CALCULATE AREA DIFFERENTIAL.
              AREA = SQRT((DQS(1)*DQT(2) - DQS(2)*DQT(1))**2
     *                  + (DQS(1)*DQT(3) - DQS(3)*DQT(1))**2
     *                  + (DQS(2)*DQT(3) - DQS(3)*DQT(2))**2)
          ELSE
C             CALCULATE ORTHOGONAL VECTOR.
              VN(1) = DQT(2)*DQS(3) - DQT(3)*DQS(2)
              VN(2) = DQT(3)*DQS(1) - DQT(1)*DQS(3)
              VN(3) = DQT(1)*DQS(2) - DQT(2)*DQS(1)
          END IF
      END IF
C
C     CALCULATE KERNEL(P,Q).
      DIST = SQRT((P(1) - Q(1))**2 + (P(2) - Q(2))**2
     *          + (P(3) - Q(3))**2)
      IF(DIST .LE. CCUTRD) THEN
C         BECAUSE OF THE NEARNESS OF P TO Q, DEFINE KERNEL(P,Q)=ZERO.
          DO 40 I=1,6
40          F(I) = cmplx(0.0,0.0)
          RETURN
      END IF
      IF(LAYER .EQ. 'SINGLE') THEN
C         CALCULATE THE SINGLE LAYER KERNEL.
          FK = AREA/DIST*(one-exp(-xk*dist))
      ELSE
C         CALCULATE THE DOUBLE LAYER KERNEL.
          DOT = VN(1)*(P(1) - Q(1)) + VN(2)*(P(2) - Q(2))
     *        + VN(3)*(P(3) - Q(3))
      FK = DOT/DIST**3*(eps*(one+xk*dist)*exp(-xk*dist)-one)
      END IF
C
C     CALCULATE INTEGRAND.
      FK = X*FK
      DO 50 I=1,6
50      F(I) = FK*C(I)
      RETURN
      END

        SUBROUTINE INTKLm(IDX,P,pn,VRINT,OPTION,LEVEL,
     *                     IFACE,VERTEX,LAYER,xk,eps)
      IMPLICIT real                       (A-H,O-Z)
      CHARACTER LAYER*6, LYR*6, OPTION*6, IOP*6
      complex xk,eps,xxk,xeps,VRINT
      DIMENSION IFACE(7,*), VERTEX(3,*), VRINT(6), VRTX(3,6), PDUP(3),
     *          P(3), VS(3), VT(3),pn(3),pnn(3)
      COMMON/PFNKLm/VRTX,AREA,PDUP,VN(3),IOP,LYR,xxk,xeps,pnn
      EXTERNAL FCNKLm
C
C     INITIALIZE COMMON FOR DEFINITION OF FCNKL2.
      IOP = OPTION
      LYR = LAYER
      xxk=xk
      xeps=eps
      DO 10 I=1,3
        pnn(i)=pn(i)
10      PDUP(I) = P(I)
      DO 20 J=1,6
        K = IFACE(J,IDX)
        DO 20 I=1,3
20        VRTX(I,J) = VERTEX(I,K)
      IF(IOP .EQ. 'PLANAR') THEN
          DO 30 I=1,3
            VT(I) = VRTX(I,2) - VRTX(I,1)
30          VS(I) = VRTX(I,3) - VRTX(I,1)
          IF(LAYER .EQ. 'SINGLE') THEN
C             CALCULATE AREA FOR PLANAR TRIANGLE.
              AREA = SQRT((VS(1)*VT(2) - VS(2)*VT(1))**2
     *                  + (VS(1)*VT(3) - VS(3)*VT(1))**2
     *                  + (VS(2)*VT(3) - VS(3)*VT(2))**2)
          ELSE
C             COMPUTE AN ORTHOGONAL VECTOR.
              VN(1) = VT(2)*VS(3) - VT(3)*VS(2)
              VN(2) = VT(3)*VS(1) - VT(1)*VS(3)
              VN(3) = VT(1)*VS(2) - VT(2)*VS(1)
          END IF
      END IF
C
C     PERFORM INTEGRATION.
      CALL SMPLX6m(VRINT,FCNKLm,LEVEL)
      RETURN
      END

      SUBROUTINE FCNKLm(S,T,F)
C     -----------------
      IMPLICIT real                       (A-H,O-Z)
      CHARACTER LAYER*6, OPTION*6
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0,
     *           FOUR=4.0D0)
      complex xk,eps,F,FK
      DIMENSION V(3,6), P(3), Q(3), F(6), C(6), DQS(3), DQT(3),pn(3)
      COMMON/PFNKLm/V,AREA,P,VN(3),OPTION,LAYER,xk,eps,pn
C
C     THE FOLLOWING CONSTANT IS TO BE A SMALL NUMBER, OF THE ORDER
C     OF THE UNIT ROUND, SAY 1000 TIMES THE UNIT ROUND.
c      CCUTRD = 1000*D1MACH(3)
        CCUTRD = 1000.0*1.81899e-32
C
C     CALCULATE CARDINAL FUNCTIONS C(I,S,T).
      SPT = S + T
      U = ONE - SPT
      C(1) = U*(TWO*U - ONE)
      C(2) = T*(TWO*T - ONE)
      C(3) = S*(TWO*S - ONE)
      C(4) = FOUR*T*U
      C(5) = FOUR*S*T
      C(6) = FOUR*S*U
C
C     CALCULATE Q=Q(S,T).
      DO 20 I=1,3
        SUM = ZERO
        DO 10 J=1,6
10        SUM = SUM + C(J)*V(I,J)
20      Q(I) = SUM
C     DO SURFACE ELEMENT AREA DIFFERENTIAL OR ORTHOGONAL VECTOR.
      IF(OPTION .EQ. 'CURVED') THEN
          DO 30 I=1,3
            DQS(I) = (FOUR*SPT - THREE)*V(I,1) + (FOUR*S - ONE)*V(I,3)
     *             + FOUR*(-T*V(I,4) + T*V(I,5)
     *             + (ONE - T - TWO*S)*V(I,6))
30          DQT(I) = (FOUR*SPT - THREE)*V(I,1) + (FOUR*T - ONE)*V(I,2)
     *             + FOUR*((ONE - S - TWO*T)*V(I,4)
     *             + S*(V(I,5) - V(I,6)))
          IF(LAYER .EQ. 'SINGLE') THEN
C             CALCULATE AREA DIFFERENTIAL.
              AREA = SQRT((DQS(1)*DQT(2) - DQS(2)*DQT(1))**2
     *                  + (DQS(1)*DQT(3) - DQS(3)*DQT(1))**2
     *                  + (DQS(2)*DQT(3) - DQS(3)*DQT(2))**2)
          ELSE
C             CALCULATE ORTHOGONAL VECTOR.
              VN(1) = DQT(2)*DQS(3) - DQT(3)*DQS(2)
              VN(2) = DQT(3)*DQS(1) - DQT(1)*DQS(3)
              VN(3) = DQT(1)*DQS(2) - DQT(2)*DQS(1)
          END IF
      END IF
C
C     CALCULATE KERNEL(P,Q).
      DIST = SQRT((P(1) - Q(1))**2 + (P(2) - Q(2))**2
     *          + (P(3) - Q(3))**2)
      IF(DIST .LE. CCUTRD) THEN
C         BECAUSE OF THE NEARNESS OF P TO Q, DEFINE KERNEL(P,Q)=ZERO.
          DO 40 I=1,6
40          F(I) = cmplx(0.0,0.0)
          RETURN
      END IF
      IF(LAYER .EQ. 'SINGLE') THEN
C         CALCULATE THE SINGLE LAYER KERNEL.
          FK = AREA/DIST*(one-exp(-xk*dist))
      ELSE
C         CALCULATE THE DOUBLE LAYER KERNEL.
          COSQ_UN = (VN(1)*(P(1) - Q(1)) + VN(2)*(P(2) - Q(2))
     *        + VN(3)*(P(3) - Q(3)))/dist
      COSP = (PN(1)*(P(1) - Q(1)) + PN(2)*(P(2) - Q(2))
     *        + PN(3)*(P(3) - Q(3)))/DIST
             xnn=vn(1)*PN(1)+vn(2)*PN(2)+vn(3)*PN(3)
       FK = (xnn-3.0*COSQ_UN*COSP)/DIST**3
      FK=FK*((one+xk*dist)*exp(-xk*dist)-one)-
     1            xk**2*exp(-xk*dist)*COSQ_UN*COSP/dist

      END IF
C
C     CALCULATE INTEGRAND.
      DO 50 I=1,6
50      F(I) = FK*C(I)
      RETURN
      END

      SUBROUTINE INTKLmm(IDXF,IDXP,pn,VRINT,OPTION,NINTEG,LAYER,
     *                  IFACE,VERTEX,xk,eps)
C
      IMPLICIT real        (A-H,O-Z)
      CHARACTER OPTION*6, IOP*6, LAYER*6, LYR*6
      PARAMETER (TWO=2.0D0)
      complex xk,eps,xxk,xeps,VRINT,VITEMP
      DIMENSION IFACE(7,*), VERTEX(3,*), VRINT(6), VRTX(3,6), P(3),
     *          VS(3), VT(3), VITEMP(6), VN(3),pn(3),pnn(3)
      COMMON/PFNKLmm/VRTX,AREA,P,pnn,IOP,LYR,IPOINT,VN,xxk,xeps
      EXTERNAL FCNKLmm
C
C     INITIALIZE COMMON FOR DEFINITION OF FCNKL4.
      IOP = OPTION
      LYR = LAYER
      xxk=xk
      xeps=eps
      DO 10 I=1,3
        pnn(i)=pn(i)
10      P(I) = VERTEX(I,IDXP)
      DO 20 J=1,6
        K = IFACE(J,IDXF)
        DO 20 I=1,3
20        VRTX(I,J) = VERTEX(I,K)
      IF(OPTION .EQ. 'PLANAR') THEN
          DO 30 I=1,3
            VT(I) = VRTX(I,2) - VRTX(I,1)
30          VS(I) = VRTX(I,3) - VRTX(I,1)
          IF(LAYER .EQ. 'SINGLE') THEN
C             CALCULATE AREA FOR PLANAR TRIANGLE.
              AREA = SQRT((VS(1)*VT(2) - VS(2)*VT(1))**2
     *                  + (VS(1)*VT(3) - VS(3)*VT(1))**2
     *                  + (VS(2)*VT(3) - VS(3)*VT(2))**2)
          ELSE
C             COMPUTE AN ORTHOGONAL VECTOR.
              VN(1) = VT(2)*VS(3) - VT(3)*VS(2)
              VN(2) = VT(3)*VS(1) - VT(1)*VS(3)
              VN(3) = VT(1)*VS(2) - VT(2)*VS(1)
          END IF
      END IF
C
C     LOCATE INDEX OF P RELATIVE TO THE FACE #IDXF.
      DO 40 I=1,6
        IF(IDXP .EQ. IFACE(I,IDXF)) THEN
            IDXP2 = I
            GO TO 45
        END IF
40      CONTINUE
      PRINT *, ' ERROR IN SUBROUTINE INTKL4.'
      PRINT *, ' THE POINT P IS NOT IN THE GIVEN TRIANGULAR FACE.'
      STOP
C
C     PERFORM INTEGRATION.
45    IF(IDXP2 .LT. 4) THEN
          IPOINT = IDXP2
          CALL INTEGRm(VRINT,FCNKLmm,NINTEG)
          RETURN
      ELSE
          IPOINT = 2*IDXP2 - 4
          CALL INTEGRm(VITEMP,FCNKLmm,NINTEG)
          IPOINT = IPOINT + 1
          CALL INTEGRm(VRINT,FCNKLmm,NINTEG)
          DO 50 J=1,6
50          VRINT(J) = (VRINT(J) + VITEMP(J))/TWO
          RETURN
      END IF
      END

      SUBROUTINE FCNKLmm(X,Y,F)
C     -----------------
      IMPLICIT real      (A-H,O-Z)
      CHARACTER OPTION*6, LAYER*6
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0,
     *           FOUR=4.0D0)
      complex xk,eps,F,FK
      DIMENSION V(3,6),P(3),Q(3),F(6),C(6),DQS(3),DQT(3),VN(3),pn(3)
      COMMON/PFNKLmm/V,AREA,P,pn,OPTION,LAYER,IPOINT,VN,xk,eps
C
C     THE FOLLOWING CONSTANT IS TO BE A SMALL NUMBER, OF THE ORDER
C     OF THE UNIT ROUND , SAY 1000 TIMES THE UNIT ROUND.
c      CCUTRD = 1000*D1MACH(3)
        CCUTRD = 1000.0*1.81899e-32
C
C     CALCULATE (S,T) FROM (X,Y).
      GO TO (1,2,3,4,5,6,7,8,9), IPOINT
1     S = Y*X
      T = (ONE - Y)*X
      GO TO 10
2     S = Y*X
      T = ONE - X
      GO TO 10
3     S = ONE - X
      T = (ONE - Y)*X
      GO TO 10
4     S = Y*X
      T = (ONE - X)/TWO + (ONE - Y)*X
      GO TO 10
5     S = (ONE - Y)*X
      T = (ONE - X)/TWO
      GO TO 10
6     S = (ONE - X)/TWO
      T = S + (ONE - Y)*X
      GO TO 10
7     T = (ONE - X)/TWO
      S = T + Y*X
      GO TO 10
8     S = (ONE - X)/TWO
      T = Y*X
      GO TO 10
9     S = (ONE - X)/TWO + Y*X
      T = (ONE - Y)*X
C
C     CALCULATE CARDINAL FUNCTIONS C(*,S,T).
10    SPT = S + T
      U = ONE - SPT
      C(1) = U*(TWO*U - ONE)
      C(2) = T*(TWO*T - ONE)
      C(3) = S*(TWO*S - ONE)
      C(4) = FOUR*T*U
      C(5) = FOUR*S*T
      C(6) = FOUR*S*U
C
C     CALCULATE Q=Q(S,T).
      DO 20 I=1,3
        SUM = ZERO
        DO 15 J=1,6
15        SUM = SUM + C(J)*V(I,J)
20      Q(I) = SUM
C
C     DO SURFACE ELEMENT AREA DIFFERENTIAL OR ORTHOGONAL VECTOR.
      IF(OPTION .EQ. 'CURVED') THEN
          DO 30 I=1,3
            DQS(I) = (FOUR*SPT - THREE)*V(I,1) + (FOUR*S - ONE)*V(I,3)
     *             + FOUR*(-T*V(I,4) + T*V(I,5)
     *             + (ONE - T - TWO*S)*V(I,6))
30          DQT(I) = (FOUR*SPT - THREE)*V(I,1) + (FOUR*T - ONE)*V(I,2)
     *             + FOUR*((ONE - S - TWO*T)*V(I,4)
     *             + S*V(I,5) - S*V(I,6))
          IF(LAYER .EQ. 'SINGLE') THEN
C             CALCULATE AREA DIFFERENTIAL.
              AREA = SQRT((DQS(1)*DQT(2) - DQS(2)*DQT(1))**2
     *                  + (DQS(1)*DQT(3) - DQS(3)*DQT(1))**2
     *                  + (DQS(2)*DQT(3) - DQS(3)*DQT(2))**2)
          ELSE
C             CALCULATE ORTHOGONAL VECTOR.
              VN(1) = DQT(2)*DQS(3) - DQT(3)*DQS(2)
              VN(2) = DQT(3)*DQS(1) - DQT(1)*DQS(3)
              VN(3) = DQT(1)*DQS(2) - DQT(2)*DQS(1)
          END IF
      END IF
C
C     CALCULATE KERNEL(P,Q).
      DIST = SQRT((P(1) - Q(1))**2 + (P(2) - Q(2))**2
     *          + (P(3) - Q(3))**2)
      IF(DIST .LE. CCUTRD) THEN
C         BECAUSE OF THE NEARNESS OF P TO Q, DEFINE KERNEL(P,Q)=ZERO.
          DO 40 I=1,6
40          F(I) = cmplx(0.0,0.0)
          RETURN
      END IF
      IF(LAYER .EQ. 'SINGLE') THEN
C         CALCULATE THE SINGLE LAYER KERNEL.
          FK = AREA/DIST*(one-exp(-xk*dist))
      ELSE
C         CALCULATE THE DOUBLE LAYER KERNEL.
          COSQ_UN = (VN(1)*(P(1) - Q(1)) + VN(2)*(P(2) - Q(2))
     *        + VN(3)*(P(3) - Q(3)))/dist
      COSP = (PN(1)*(P(1) - Q(1)) + PN(2)*(P(2) - Q(2))
     *        + PN(3)*(P(3) - Q(3)))/DIST
             xnn=vn(1)*PN(1)+vn(2)*PN(2)+vn(3)*PN(3)
       FK = (xnn-3.0*COSQ_UN*COSP)/DIST**3
      FK=FK*((one+xk*dist)*exp(-xk*dist)-one)-
     1            xk**2*exp(-xk*dist)*COSQ_UN*COSP/dist
      END IF
C
C     CALCULATE INTEGRAND.
      FK = X*FK
      DO 50 I=1,6
50      F(I) = FK*C(I)
      RETURN
      END

      SUBROUTINE INTKL1m(IDX,P,VRINT,OPTION,LEVEL,IFACE,VERTEX)
      IMPLICIT real (A-H,O-Z)
      CHARACTER OPTION*6,IOP*6
      complex VRINT
      DIMENSION IFACE(7,*), VERTEX(3,*), VRINT(6), VRTX(3,6), PDUP(3),
     *          P(3), VS(3), VT(3)
      COMMON/PFNKL1/VRTX,AREA,PDUP,IOP
      EXTERNAL FCNKL1
C
C     INITIALIZE COMMON FOR DEFINITION OF FCNKL1.
      IOP = OPTION
      DO 10 I=1,3
10      PDUP(I) = P(I)
      DO 20 J=1,6
        K = IFACE(J,IDX)
        DO 20 I=1,3
20        VRTX(I,J) = VERTEX(I,K)
      IF(OPTION .EQ. 'PLANAR') THEN
          DO 30 I=1,3
            VT(I) = VRTX(I,2) - VRTX(I,1)
30          VS(I) = VRTX(I,3) - VRTX(I,1)
C         CALCULATE AREA OF PLANAR TRIANGLE.
          AREA = SQRT((VS(1)*VT(2) - VS(2)*VT(1))**2
     *              + (VS(1)*VT(3) - VS(3)*VT(1))**2
     *              + (VS(2)*VT(3) - VS(3)*VT(2))**2)
      END IF
C
C     PERFORM INTEGRATION.
      CALL SMPLX6m(VRINT,FCNKL1,LEVEL)
      RETURN
      END

        SUBROUTINE FCNKL1(S,T,F)
      IMPLICIT real        (A-H,O-Z)
      complex        KERNEL,F,FK
      CHARACTER OPTION*6
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0,
     *           FOUR=4.0D0)
      DIMENSION V(3,6), P(3), Q(3), F(6), C(6), DQS(3), DQT(3)
      COMMON/PFNKL1/V,AREA,P,OPTION
C
C     CALCULATE CARDINAL FUNCTIONS C(*,S,T).
      SPT = S + T
      U = ONE - SPT
      C(1) = U*(TWO*U - ONE)
      C(2) = T*(TWO*T - ONE)
      C(3) = S*(TWO*S - ONE)
      C(4) = FOUR*T*U
      C(5) = FOUR*S*T
      C(6) = FOUR*S*U
C
C     CALCULATE Q=Q(S,T).
      DO 20 I=1,3
        SUM = ZERO
        DO 10 J=1,6
10        SUM = SUM + C(J)*V(I,J)
20      Q(I) = SUM
C
C     DO SURFACE ELEMENT AREA DIFFERENTIAL.
      IF(OPTION .EQ. 'CURVED') THEN
          DO 30 I=1,3
            DQS(I) = (FOUR*SPT - THREE)*V(I,1) + (FOUR*S - ONE)*V(I,3)
     *             + FOUR*(-T*V(I,4) + T*V(I,5)
     *             + (ONE - T - TWO*S)*V(I,6))
30          DQT(I) = (FOUR*SPT - THREE)*V(I,1) + (FOUR*T - ONE)*V(I,2)
     *             + FOUR*((ONE - S - TWO*T)*V(I,4)
     *             + S*(V(I,5) - V(I,6)))
C         CALCULATE AREA DIFFERENTIAL.
          AREA = SQRT((DQS(1)*DQT(2) - DQS(2)*DQT(1))**2
     *              + (DQS(1)*DQT(3) - DQS(3)*DQT(1))**2
     *              + (DQS(2)*DQT(3) - DQS(3)*DQT(2))**2)
      END IF
C
C     CALCULATE INTEGRAND.
      FK = KERNEL(P,Q)*AREA
      DO 50 I=1,6
50      F(I) = FK*C(I)
      RETURN
      END

      SUBROUTINE INTKL3m(IDXF,IDXP,VRINT,OPTION,NINTEG,IFACE,VERTEX)
      IMPLICIT real   (A-H,O-Z)
      CHARACTER OPTION*6, IOP*6
      PARAMETER (TWO=2.00)
      complex  VRINT,VITEMP
      DIMENSION IFACE(7,*), VERTEX(3,*), VRINT(6), VRTX(3,6), P(3),
     *          VS(3), VT(3), VITEMP(6)
      COMMON/PFNKL3/VRTX,AREA,P,IPOINT,IOP
      EXTERNAL FCNKL3
C
C     INITIALIZE COMMON FOR DEFINITION OF FCNKL3.
      IOP = OPTION
      DO 10 I=1,3
10      P(I) = VERTEX(I,IDXP)
      DO 20 J=1,6
        K = IFACE(J,IDXF)
        DO 20 I=1,3
20        VRTX(I,J) = VERTEX(I,K)
      IF(OPTION .EQ. 'PLANAR') THEN
          DO 30 I=1,3
            VT(I) = VRTX(I,2) - VRTX(I,1)
30          VS(I) = VRTX(I,3) - VRTX(I,1)
C         CALCULATE AREA OF PLANAR TRIANGLE.
          AREA = SQRT((VS(1)*VT(2) - VS(2)*VT(1))**2
     *              + (VS(1)*VT(3) - VS(3)*VT(1))**2
     *              + (VS(2)*VT(3) - VS(3)*VT(2))**2)
      END IF
C
C     LOCATE INDEX OF P RELATIVE TO THE FACE #IDXF.
      DO 40 I=1,6
        IF(IDXP .EQ. IFACE(I,IDXF)) THEN
            IDXP2 = I
            GO TO 45
        END IF
40      CONTINUE
      PRINT *, ' ERROR IN SUBROUTINE INTKL3.'
      PRINT *, ' THE POINT P IS NOT IN THE GIVEN TRIANGULAR FACE.'
      STOP
C
C     PERFORM INTEGRATION.
45    IF(IDXP2 .LT. 4) THEN
          IPOINT = IDXP2
          CALL INTEGRm(VRINT,FCNKL3,NINTEG)
          RETURN
      ELSE
          IPOINT = 2*IDXP2 - 4
          CALL INTEGRm(VITEMP,FCNKL3,NINTEG)
          IPOINT = IPOINT + 1
          CALL INTEGRm(VRINT,FCNKL3,NINTEG)
          DO 50 J=1,6
50          VRINT(J) = (VRINT(J) + VITEMP(J))/TWO
          RETURN
      END IF
      END

      SUBROUTINE FCNKL3(X,Y,F)
      IMPLICIT real      (A-H,O-Z)
      complex       KERNEL,F,FK
      CHARACTER OPTION*6
      PARAMETER (ZERO=0.00, ONE=1.00, TWO=2.00, THREE=3.00,
     *           FOUR=4.00)
      DIMENSION V(3,6), P(3), Q(3), F(6), C(6), DQS(3), DQT(3)
      COMMON/PFNKL3/V,AREA,P,IPOINT,OPTION
C
C     CALCULATE (S,T) FROM (X,Y).
      GO TO (1,2,3,4,5,6,7,8,9), IPOINT
1     S = Y*X
      T = (ONE - Y)*X
      GO TO 10
2     S = Y*X
      T = ONE - X
      GO TO 10
3     S = ONE - X
      T = (ONE - Y)*X
      GO TO 10
4     S = Y*X
      T = (ONE - X)/TWO + (ONE - Y)*X
      GO TO 10
5     S = (ONE - Y)*X
      T = (ONE - X)/TWO
      GO TO 10
6     S = (ONE - X)/TWO
      T = S + (ONE - Y)*X
      GO TO 10
7     T = (ONE - X)/TWO
      S = T + Y*X
      GO TO 10
8     S = (ONE - X)/TWO
      T = Y*X
      GO TO 10
9     S = (ONE - X)/TWO + Y*X
      T = (ONE - Y)*X
C
C     CALCULATE CARDINAL FUNCTIONS C(*,S,T).
10    SPT = S + T
      U = ONE - SPT
      C(1) = U*(TWO*U - ONE)
      C(2) = T*(TWO*T - ONE)
      C(3) = S*(TWO*S - ONE)
      C(4) = FOUR*T*U
      C(5) = FOUR*S*T
      C(6) = FOUR*S*U
C
C     CALCULATE Q=Q(S,T).
      DO 20 I=1,3
        SUM = ZERO
        DO 15 J=1,6
15        SUM = SUM + C(J)*V(I,J)
20      Q(I) = SUM
C
C     DO SURFACE ELEMENT AREA DIFFERENTIAL.
      IF(OPTION .EQ. 'CURVED') THEN
          DO 30 I=1,3
            DQS(I) = (FOUR*SPT - THREE)*V(I,1) + (FOUR*S - ONE)*V(I,3)
     *             + FOUR*(-T*V(I,4) + T*V(I,5)
     *             + (ONE - T - TWO*S)*V(I,6))
30          DQT(I) = (FOUR*SPT - THREE)*V(I,1) + (FOUR*T - ONE)*V(I,2)
     *             + FOUR*((ONE - S - TWO*T)*V(I,4)
     *             + S*V(I,5) - S*V(I,6))
C         CALCULATE AREA DIFFERENTIAL.
          AREA = SQRT((DQS(1)*DQT(2) - DQS(2)*DQT(1))**2
     *              + (DQS(1)*DQT(3) - DQS(3)*DQT(1))**2
     *              + (DQS(2)*DQT(3) - DQS(3)*DQT(2))**2)
      END IF
C
C     CALCULATE INTEGRAND.
      FK = X*KERNEL(P,Q)*AREA
      DO 50 I=1,6
50      F(I) = FK*C(I)
      RETURN
      END

           SUBROUTINE SMPLX6m(VRINT,FCNKL,LEVEL)
      IMPLICIT none
      real zero,one,two,sr15,a,b,c,t,r,s,u,v
      PARAMETER (ZERO=0.0, ONE=1.0, TWO=2.0,
     *           SR15=3.87298334620741)
      PARAMETER (A=.1125, B=(155.0-SR15)/2400.0,
     *           C=(155.0+SR15)/2400.0, T=ONE/3.0,
     *           R=(6.0-SR15)/21.0, S=(9.0+TWO*SR15)/21.0,
     *           U=(6.0+SR15)/21.0, V=(9.0-TWO*SR15)/21.0)
       complex VRINT,FN
      real sm,tm,w,sv,tv,sign,h
      integer i,j,k,l,nr,nf,level
      DIMENSION SM(7),TM(7),W(7),FN(6),VRINT(6)
C
C     INITIALIZE.
      SM(1) = T
      TM(1) = T
      SM(2) = R
      TM(2) = R
      SM(3) = S
      TM(3) = R
      SM(4) = R
      TM(4) = S
      SM(5) = U
      TM(5) = U
      SM(6) = V
      TM(6) = U
      SM(7) = U
      TM(7) = V
      W(1) = A
      DO 10 I=1,3
        W(I+1) = B
10      W(I+4) = C
C
C     INITIALIZE VARIABLES FOR INTEGRATION LOOP.
      NF = 4**LEVEL
      NR = 2**LEVEL
      H = ONE/NR
      IF(LEVEL .GT. 0) THEN
          DO 20 K=1,7
            SM(K) = SM(K)/NR
20          TM(K) = TM(K)/NR
      END IF
      DO 25 L=1,6
25      VRINT(L) = cmplx(0.0,0.0)
C
C     BEGIN INTEGRATION OVER NF SUB-SIMPLICES OF UNIT SIMPLEX.
      DO 40 I=1,NR
        SV = ZERO
        TV =(NR - I)*H
        SIGN = ONE
        DO 40 J=1,2*I-1
          DO 30 K=1,7
            CALL FCNKL(SV+SIGN*SM(K), TV+SIGN*TM(K), FN)
            DO 30 L=1,6
30            VRINT(L) = VRINT(L) + W(K)*FN(L)
          SIGN = -SIGN
          IF(SIGN .LT. ZERO) THEN
              SV = SV + H
              TV = TV + H
          ELSE
              TV = TV - H
          END IF
40        CONTINUE
      DO 50 L=1,6
50      VRINT(L) = VRINT(L)/NF
      RETURN
      END

        SUBROUTINE INTEGRm(VINT,F,N)
      IMPLICIT none
      real zero,one
      PARAMETER (ZERO=0.0, ONE=1.0)
      complex VINT,FVEC,SUM
       real w,t
       integer i,j,k,n
      DIMENSION W(20), T(20), VINT(6), SUM(6), FVEC(6)
C
      CALL ZEROLG(N,T,W,ZERO,ONE)
      DO 5 K=1,6
5        VINT(K) = cmplx(0.0,0.0)
      DO 30 I=1,N
        DO 10 K=1,6
10        SUM(K) = cmplx(0.0,0.0)
        DO 20 J=1,N
          CALL F(T(I),T(J),FVEC)
          DO 20 K=1,6
20          SUM(K) = SUM(K) + W(J)*FVEC(K)
        DO 30 K=1,6
30        VINT(K) = VINT(K) + W(I)*SUM(K)
      RETURN
      END

      SUBROUTINE ZEROLG(N,ZZ,WW,A,B) 
      IMPLICIT real  (A-H,O-Z)
      DIMENSION Z(109), W(109), ZZ(N), WW(N)
      DATA ONE/1.0/, TWO/2.0/
      DATA Z(1),Z(2),Z(3),Z(4),Z(5),Z(6),Z(7),Z(8),Z(9),Z(10)/
     *.577350269189626D0,.774596669241483D0,0.0D0,.861136311594053D0,
     *.339981043584856D0,.906179845938664D0,.538469310105683D0,0.0D0,
     *.932469514203152D0,.661209386466265D0/
      DATA Z(11),Z(12),Z(13),Z(14),Z(15),Z(16),Z(17),Z(18),Z(19),Z(20)/
     *.238619186083197D0,.949107912342759D0,.741531185599394D0,
     *.405845151377397D0,0.0D0,.960289856497536D0,.796666477413627D0,
     *.525532409916329D0,.183434642495650D0,.968160239507626D0/
      DATA Z(21),Z(22),Z(23),Z(24),Z(25),Z(26),Z(27),Z(28),Z(29)/
     *.836031107326636D0,.613371432700590D0,.324253423403809D0,
     *0.0D0,.973906528517172D0,.865063366688985D0,.679409568299024D0,
     *.433395394129247D0,.148874338981631D0/
      DATA Z(30),Z(31),Z(32),Z(33),Z(34),Z(35),Z(36),Z(37),Z(38)/
     *.978228658146057D0,.887062599768095D0,.730152005574049D0,
     *.519096129206812D0,.269543155952345D0,0.0D0,
     *.981560634246719D0,.904117256370475D0,.769902674194305D0/
      DATA Z(39),Z(40),Z(41),Z(42),Z(43),Z(44),Z(45),Z(46),Z(47)/
     *.587317954286617D0,.367831498998180D0,.125233408511469D0,
     *.984183054718588D0,.917598399222978D0,.801578090733310D0,
     *.642349339440340D0,.448492751036447D0,.230458315955135D0/
      DATA Z(48),Z(49),Z(50),Z(51),Z(52),Z(53),Z(54),Z(55),Z(56)/
     *0.0D0,.986283808696812D0,.928434883663574D0,
     *.827201315069765D0,.687292904811685D0,.515248636358154D0,
     *.319112368927890D0,.108054948707344D0,.987992518020485D0/
      DATA Z(57),Z(58),Z(59),Z(60),Z(61),Z(62),Z(63),Z(64),Z(65)/
     *.937273392400706D0,.848206583410427D0,.724417731360170D0,
     *.570972172608539D0,.394151347077563D0,.201194093997435D0,
     *0.0D0,.989400934991650D0,.944575023073233D0/
      DATA Z(66),Z(67),Z(68),Z(69),Z(70),Z(71),Z(72),Z(73),Z(74)/
     *.865631202387832D0,.755404408355003D0,.617876244402644D0,
     *.458016777657227D0,.281603550779259D0,.0950125098376374D0,
     *.990575475314417D0,.950675521768768D0,.880239153726986D0/
      DATA Z(75),Z(76),Z(77),Z(78),Z(79),Z(80),Z(81),Z(82),Z(83)/
     *.781514003896801D0,.657671159216691D0,.512690537086477D0,
     *.351231763453876D0,.178484181495848D0,0.0D0,
     *.991565168420931D0,.955823949571398D0,.892602466497556D0/
      DATA Z(84),Z(85),Z(86),Z(87),Z(88),Z(89),Z(90),Z(91),Z(92)/
     *.803704958972523D0,.691687043060353D0,.559770831073948D0,
     *.411751161462843D0,.251886225691506D0,.0847750130417353D0,
     *.992406843843584D0,.960208152134830D0,.903155903614818D0/
      DATA Z(93),Z(94),Z(95),Z(96),Z(97),Z(98),Z(99),Z(100),Z(101)/
     *.822714656537143D0,.720966177335229D0,.600545304661681D0,
     *.464570741375961D0,.316564099963630D0,.160358645640225D0,
     *0.0D0,.993128599185095D0,.963971927277914D0/
      DATA Z(102),Z(103),Z(104),Z(105),Z(106),Z(107),Z(108),Z(109)/
     *.912234428251326D0,.839116971822219D0,.746331906460151D0,
     *.636053680726515D0,.510867001950827D0,.373706088715420D0,
     *.227785851141645D0,.0765265211334973D0/
      DATA W(1),W(2),W(3),W(4),W(5),W(6),W(7),W(8),W(9),W(10)/
     *1.0D0,.555555555555556D0,.888888888888889D0,.347854845137454D0,
     *.652145154862546D0,.236926885056189D0,.478628670499366D0,
     *.568888888888889D0,.171324492379170D0,.360761573048139D0/
      DATA W(11),W(12),W(13),W(14),W(15),W(16),W(17),W(18),W(19),W(20)/
     *.467913934572691D0,.129484966168870D0,.279705391489277D0,
     *.381830050505119D0,.417959183673469D0,.101228536290376D0,
     *.222381034453374D0,.313706645877887D0,.362683783378362D0,
     *.0812743883615744D0/
      DATA W(21),W(22),W(23),W(24),W(25),W(26),W(27),W(28),W(29)/
     *.180648160694857D0,.260610696402935D0,.312347077040003D0,
     *.330239355001260D0,.0666713443086881D0,.149451349150581D0,
     *.219086362515982D0,.269266719309996D0,.295524224714753D0/
      DATA W(30),W(31),W(32),W(33),W(34),W(35),W(36),W(37),W(38)/
     *.0556685671161737D0,.125580369464905D0,.186290210927734D0,
     *.233193764591990D0,.262804544510247D0,.272925086777901D0,
     *.0471753363865118D0,.106939325995318D0,.160078328543346D0/
      DATA W(39),W(40),W(41),W(42),W(43),W(44),W(45),W(46),W(47)/
     *.203167426723066D0,.233492536538355D0,.249147045813403D0,
     *.0404840047653159D0,.0921214998377284D0,.138873510219787D0,
     *.178145980761946D0,.207816047536889D0,.226283180262897D0/
      DATA W(48),W(49),W(50),W(51),W(52),W(53),W(54),W(55),W(56)/
     *.232551553230874D0,.0351194603317519D0,.0801580871597602D0,
     *.121518570687903D0,.157203167158194D0,.185538397477938D0,
     *.205198463721296D0,.215263853463158D0,.0307532419961173D0/
      DATA W(57),W(58),W(59),W(60),W(61),W(62),W(63),W(64),W(65)/
     *.0703660474881081D0,.107159220467172D0,.139570677926154D0,
     *.166269205816994D0,.186161000015562D0,.198431485327112D0,
     *.202578241925561D0,.0271524594117541D0,.0622535239386478D0/
      DATA W(66),W(67),W(68),W(69),W(70),W(71),W(72),W(73),W(74)/
     *.0951585116824928D0,.124628971255534D0,.149595988816577D0,
     *.169156519395003D0,.182603415044924D0,.189450610455068D0,
     *.0241483028685479D0,.0554595293739872D0,.0850361483171792D0/
      DATA W(75),W(76),W(77),W(78),W(79),W(80),W(81),W(82),W(83)/
     *.111883847193404D0,.135136368468525D0,.154045761076810D0,
     *.168004102156450D0,.176562705366993D0,.179446470356207D0,
     *.0216160135264833D0,.0497145488949698D0,.0764257302548891D0/
      DATA W(84),W(85),W(86),W(87),W(88),W(89),W(90),W(91),W(92)/
     *.100942044106287D0,.122555206711478D0,.140642914670651D0,
     *.154684675126265D0,.164276483745833D0,.169142382963144D0,
     *.0194617882297265D0,.0448142267656996D0,.0690445427376412D0/
      DATA W(93),W(94),W(95),W(96),W(97),W(98),W(99),W(100),W(101)/
     *.0914900216224500D0,.111566645547334D0,.128753962539336D0,
     *.142606702173607D0,.152766042065860D0,.158968843393954D0,
     *.161054449848784D0,.0176140071391521D0,.0406014298003869D0/
      DATA W(102),W(103),W(104),W(105),W(106),W(107),W(108),W(109)/
     *.0626720483341091D0,.0832767415767047D0,.101930119817240D0,
     *.118194531961518D0,.131688638449177D0,.142096109318382D0,
     *.149172986472604D0,.152753387130726D0/
C
C     CALCULATE THE WEIGHTS AND NODES.
C
      IF(N .EQ. 1) THEN
          ZZ(1) = (A + B)/TWO
          WW(1) = B - A
          RETURN
      END IF
C
      SCALE = (B - A)/TWO
      IF((N/2)*2 .NE. N) THEN
          IBASE = (N*N - 1)/4
          M = (N - 1)/2
          ZZ(M+1) = (A + B)/TWO
          WW(M+1) = W(IBASE+M)*SCALE
      ELSE
          M = N/2
          IBASE = M*M
      END IF
C
      DO 120 I=1,M
        T = Z(IBASE+I-1)
        ZZ(I) = (A*(ONE + T) + (ONE - T)*B)/TWO
        ZZ(N+1-I) = (A*(ONE - T) + (ONE + T)*B)/TWO
        WW(I) = W(IBASE+I-1)*SCALE
120     WW(N+1-I) = WW(I)
      RETURN
      END


