c$$$! All variables made lower case (JH) 4/26/05
c$$$
c$$$C*********************************************************************C
c$$$C*                                                                   *C
c$$$C*  chapman.for                                                      *C
c$$$C*                                                                   *C
c$$$C*  Written by:  David L. Huestis, Molecular Physics Laboratory      *C
c$$$C*                                                                   *C
c$$$C*  Copyright (c) 2000  SRI International                            *C
c$$$C*  All Rights Reserved                                              *C
c$$$C*                                                                   *C
c$$$C*  This software is provided on an as is basis; without any         *C
c$$$C*  warranty; without the implied warranty of merchantability or     *C
c$$$C*  fitness for a particular purpose.                                *C
c$$$C*                                                                   *C
c$$$C*********************************************************************C
c$$$C*
c$$$C*	To calculate the Chapman Function, Ch(X,chi0), the column 
c$$$C*	depth of an exponential atmosphere integrated along a line 
c$$$C*	from a given point to the sun, divided by the column depth for 
c$$$C*	a vertical sun.
c$$$C*
c$$$C*  USAGE:
c$$$C*
c$$$C*	  z = altitude above the surface
c$$$C*	  R = radius of the planet
c$$$C*	  H = atmospheric scale height
c$$$C*
c$$$C*	  X = (R+z)/H
c$$$C*	  chi0 = solar zenith angle (in degrees)
c$$$C*
c$$$C*	  implicit real*4(a-h,o-z)
c$$$C*	  depth = atm_chapman(X,chi0)	! analytical
c$$$C*	  depth = atm_chap_num(X,chi0)	! numerical (chi0 .le. 90)
c$$$C*
c$$$C*	  implicit real*8(a-h,o-z)
c$$$C*	  depth = atm8_chapman(X,chi0)	! analytical
c$$$C*	  depth = atm8_chap_num(X,chi0)	! numerical (chi0 .le. 90)
c$$$C*
c$$$C*  PERFORMANCE:
c$$$C*
c$$$C*	Compiled and linked using Microsoft FORTRAN 5.1, and executed 
c$$$C*	in MS-DOS mode under Windows 95 on a 160 MHz PC.
c$$$C*
c$$$C*    TIMING (in microseconds, typical)
c$$$C*
c$$$C*	  120	atm_chapman and atm8_chapman for X .lt. 36
c$$$C*	   25	atm_chapman and atm8_chapman for X .ge. 36
c$$$C*	  500	atm_chap_num
c$$$C*	 5000	atm8_chap_num
c$$$C*
c$$$C*    ACCURACY (maximum relative error, 0.le.chi0.le.90, 1.le.X.le.820)
c$$$C*
c$$$C*	6.0E-7	atm_chapman and atm8_chapman for X .lt. 60
c$$$C*	1.5E-7	atm_chapman and atm8_chapman for X .ge. 60
c$$$C*	6.0E-8	atm_chap_num
c$$$C*	1.E-15	atm8_chap_num (convergence test)
c$$$C*
c$$$C*    CODING
c$$$C*
c$$$C*	No claims are made that the code is optimized for speed, 
c$$$C*	accuracy, or compactness.  The principal objectives were 
c$$$C*
c$$$C*	  (1) Robustness with respect to argument values
c$$$C*	  (2) Rigorous mathematical derivation and error control
c$$$C*	  (3) Maximal use of "well known" mathematical functions
c$$$C*	  (4) Ease of readability and mapping of theory to coding
c$$$C*
c$$$C*	The real*8 accuracy could be improved with more accurate 
c$$$C*	representations of E1(), erfc(), I0(), I1(), K0(), K1().
c$$$C*
c$$$C*	In the course of development, many representations and 
c$$$C*	approximations of the Chapman Function were attempted that 
c$$$C*	failed to be robustly extendable to machine-precision.
c$$$C*
c$$$C*  INTERNET ACCESS:
c$$$C*
c$$$C*	Source: http://www-mpl.sri.com/software/chapman/chapman.html
c$$$C*	Author: mailto:david.huestis@sri.com
c$$$C*	        http://www-mpl.sri.com/bios/Huestis-DL.html
c$$$C*
c$$$C*  EDIT HISTORY:
c$$$C*
c$$$C*	01/22/2000 DLH	First complete documentation
c$$$C*
c$$$C*	01/15/2000 DLH	First complete version of chapman.for
c$$$C*
c$$$C**********************************************************************
c$$$C*
c$$$C*  THEORY:
c$$$C*
c$$$C*    INTRODUCTION
c$$$C*
c$$$C*	    This computer code models the absorption of solar radiation 
c$$$C*	by an atmosphere that depends exponentionally on altitude.  In 
c$$$C*	specific we calculate the effective column depth of a species 
c$$$C*	of local density, n(z), from a point at a given altitude, z0, 
c$$$C*	to the sun at a given solar zenith angle, chi0.  Following Rees 
c$$$C*	[Re89, Section 2.2] we write the column depth for chi0 .le. 90 
c$$$C*	degrees as
c$$$C*
c$$$C*   (A)  N(z0,chi0) = int{z=z0,infinity} 
c$$$C*	     [ n(z)/sqrt( 1 - ( sin(chi0) * (R+z0) / (R+z) ) **2 ) dz ]
c$$$C*
c$$$C*	where R is the radius of the solid planet (e.g. Earth).  For 
c$$$C*	chi0 .gt. 90 degrees we write
c$$$C*
c$$$C*	  N(z0,chi0) = 2*N(zs,90) - N(z0,180-chi0)
c$$$C*
c$$$C*	where zs = (R+z0)*sin(chi0)-R is the tangent height.
c$$$C*
c$$$C*	    For an exponential atmosphere, with
c$$$C*
c$$$C*	  n(z) = n(z0) * exp(-(z-z0)/H)
c$$$C*
c$$$C*	with a constant scale height, H, the column depth can be 
c$$$C*	represented by the Chapman function, Ch(X,chi0), named after 
c$$$C*	the author of the first quantitative mathematical investigation 
c$$$C*	[Ch31b] trough the relation
c$$$C*
c$$$C*	  N(z0,chi0) = H * n(z0) * Ch(X,chi0)
c$$$C*
c$$$C*	where X = (R+z0)/H is a dimensionless measure of the radius 
c$$$C*	of curvature, with values from about 300 to 1300 on Earth.
c$$$C*
c$$$C*
c$$$C*    APPROACH
c$$$C*
c$$$C*	    We provide function entry points for very stable and 
c$$$C*	reasonably efficient evaluation of Ch(X,chi0) with full 
c$$$C*	single-precision accuracy (.le. 6.0E-7 relative) for a wide 
c$$$C*	range of parameters.  A 15-digit-accurate double precision 
c$$$C*	numerical integration routine is also provided.
c$$$C*
c$$$C*	    Below we will develop (1) a compact asymptotic expansion of 
c$$$C*	good accuracy for moderately large values of X (.gt. 36) and all 
c$$$C*	values of chi0, (2) an efficient numerical integral for 
c$$$C*	all values of X and chi0, and (3) an explicit analytical 
c$$$C*	representation, valid for all values of X and chi0, based 
c$$$C*	the differential equation satisfied by Ch(X,chi0).
c$$$C*
c$$$C*	    All three of these represent new research results as well 
c$$$C*	as significant computational improvements over the previous 
c$$$C*	litearture, much of which is cited below.
c$$$C*
c$$$C*
c$$$C*    CHANGES OF THE VARIABLE OF INTEGRATION
c$$$C*
c$$$C*	Substituting y = (R+z)/(R+z0) - 1 we find
c$$$C*
c$$$C*   (B)  Ch(X,chi0) = X * int{y=0,infinity}
c$$$C*	     [ exp(-X*y) / sqrt( 1 - ( sin(chi0) / (1+y) )**2 ) dy ]
c$$$C*
c$$$C*	The futher substitutions s = (1+y)/sin(chi0), s0 = 1/sin(chi0) 
c$$$C*	give
c$$$C*
c$$$C*   (C)  Ch(X,chi0) = X*sin(chi0) * int{s=s0,infinity}
c$$$C*	     [ exp(X*(1-sin(chi0)*s)) * s / sqrt(s**2-1) ds ]
c$$$C*
c$$$C*	From this equation we can establish that
c$$$C*
c$$$C*	  Ch(X,90) = X*exp(X)*K1(X)
c$$$C*
c$$$C*	[AS64, Equations 9.6.23 and 9.6.27].  If we now substitute
c$$$C*	s = 1/sin(lambda) we obtain
c$$$C*
c$$$C*   (D)  Ch(X,chi0) = X*sin(chi0) * int{lambda=0,chi0} 
c$$$C*	    [ exp(X*(1-sin(chi0)*csc(lambda))) * csc(lambda)**2 dlambda]
c$$$C*
c$$$C*	which is the same as Chapman's original formulation [Ch31b, p486,
c$$$C*	eqn (10)].  If we first expand the square root in (B)
c$$$C*
c$$$C*	  1/sqrt(1-q) = 1 + q/( sqrt(1-q)*(1+sqrt(1-q)) )
c$$$C*
c$$$C*	with q = ( sin(chi0) / (1+y) )**2 = sin(lambda)**2, we obtain 
c$$$C*	a new form of (D) without numerical sigularities and simple 
c$$$C*	convergence to Ch(0,chi0) = Ch(X,0) = 1
c$$$C*
c$$$C*   (E)  Ch(X,chi0) = 1 + X*sin(chi0) * int{lambda=0,chi0} 
c$$$C*	    [ exp(X*(1-sin(chi0)*csc(lambda))) 
c$$$C*		/ (1 + cos(lambda) ) dlambda ]
c$$$C*
c$$$C*	Alternatively, we may substitute t**2 = y + t0**2, 
c$$$C*	into Equation (B), with t0**2 = 1-sin(chi0), finding
c$$$C*
c$$$C*   (F)  Ch(X,chi0) = X * int{s=t0,infinity} 
c$$$C*	    [ exp(-X*(t**2-t0**2)) * f(t,chi0) dt ]
c$$$C* 
c$$$C*	where
c$$$C*
c$$$C*	  f(t,chi0) = (t**2 + sin(chi0)) / sqrt(t**2+2*sin(chi0))
c$$$C*
c$$$C*	  f(t,chi0) = (t**2-t0**2+1)/sqrt(t**2-t0**2+1+sin(chi0))
c$$$C*
c$$$C*	    Below we will use Equation (F) above to develop a
c$$$C*	compact asymptotic expansion of good accuracy for moderately 
c$$$C*	large values of X (.gt. 36) and all values of chi0, Equation (E) 
c$$$C*	to develop an efficient numerical integral for Ch(X,chi0) for 
c$$$C*	all values of X and chi0, and Equation (C) to derive an explicit 
c$$$C*	analytical representation, valid for all values of X and chi0,  
c$$$C*	based on the differential equation satisfied by Ch(X,chi0).
c$$$C*
c$$$C*    atm_chapman(X,chi0) and atm8_chapman(X,chi0)
c$$$C*
c$$$C*	These routines return real*4 and real*8 values of Ch(X,chi0)
c$$$C*	selecting the asymptotic expansion or differential equation 
c$$$C*	approaches, depending on the value of X.  These routines also 
c$$$C*	handle the case of chi0 .gt. 90 degrees.
c$$$C*
c$$$C*    atm_chap_num(X,chi0) and atm8_chap_num(X,chi0)
c$$$C*
c$$$C*	These routines return real*4 and real*8 values of Ch(X,chi0)
c$$$C*	evaluated numerically.  They are both more accurate than the 
c$$$C*	corresponding atm*_chapman() functions, but take significantly 
c$$$C*	more CPU time.
c$$$C*
c$$$C*
c$$$C*    ASYMPTOTIC EXPANSION
c$$$C*
c$$$C*	From Equation (F) we expand, with t0**2 = 1-sin(chi0), 
c$$$C*
c$$$C*	  f(t,chi0) = sum{n=0,3} [ C(n,chi0) * (t**2-t0**2)**n ]
c$$$C*
c$$$C*	The function atm8_chap_asy(X,chi0) evaluates integrals of the 
c$$$C*	form
c$$$C*
c$$$C*	  int{t=t0,infinity} [exp(-X*(t**2-t0**2))*(t**2-t0**2)**n dt]
c$$$C*
c$$$C*	in terms of incomplete gamma functions, and sums them to 
c$$$C*	compute Ch(X,chi0).  For large values of X, this results in an 
c$$$C*	asymptotic expansion in negative powers of X, with coefficients 
c$$$C*	that are stable for all values of chi0.
c$$$C*
c$$$C*	In contrast, the asymptotic expansions of Chapman [Ch31b, 
c$$$C*	p488, Equation (22) and p490, Equation (38)], Hulburt [He39], 
c$$$C*	and Swider [Sw64, p777, Equation (43)] use negative powers of 
c$$$C*	X*cos(chi0)**2 or X*sin(chi0), and are accurate only for 
c$$$C*	small values or large values of chi0, respectively.
c$$$C*
c$$$C*	Taking only the first term in the present expansion gives the 
c$$$C*	simple formula
c$$$C*
c$$$C*	  Ch(X,chi0) = sqrt(pi*X/(1+sin(chi0))) * exp(X*(1-sin(chi0)))
c$$$C*		* erfc( sqrt(X*(1-sin(chi0))) )
c$$$C*
c$$$C*	This is slightly more accurate than the semiempirical 
c$$$C*	formula of Fitzmaurice [Fi64, Equation (3)], and sightly less 
c$$$C*	accurate than that of Swider [Sw64, p780, Equation (52), 
c$$$C*	corrected in SG69].
c$$$C*
c$$$C*
c$$$C*    NUMERICAL INTEGRATION
c$$$C*
c$$$C*	We are integrating
c$$$C*
c$$$C*   (E)  Ch(X,chi0) = 1 + X*sin(chi0) * int{lambda=0,chi0} 
c$$$C*	    [ exp(X*(1-sin(chi0)*csc(lambda))) 
c$$$C*		/ ( 1 + cos(lambda) ) dlambda ]
c$$$C*
c$$$C*	The integrand is numerically very smooth, and rapidly varying 
c$$$C*	only near lambda = 0.  For X .ne. 0 we choose the lower limit 
c$$$C*	of numerical integration such that the integrand is 
c$$$C*	exponentially small, 7.0E-13 (3.0E-20 for real*8).  The domain 
c$$$C*	of integration is divided into 64 equal intervals (6000 for 
c$$$C*	real*8), and integrated numerically using the 9-point closed 
c$$$C*	Newton-Cotes formula from Hildebrand [Hi56a, page 75, Equation
c$$$C*	(3.5.17)].
c$$$C*
c$$$C*
c$$$C*    INHOMOGENOUS DIFFERENTIAL EQUATION
c$$$C*
c$$$C*	    The function atm8_chap_deq(X,chi0) calculates Ch(X,chi0), 
c$$$C*	based on Equation (C) above, using the inhomogeneous 
c$$$C*	Bessel's equation as described below.  Consider the function 
c$$$C*
c$$$C*	  Z(Q) = int{s=s0,infinity} [ exp(-Q*s) / sqrt(s**2-1) ds ]
c$$$C*
c$$$C*	Differentiating with respect to Q we find that 
c$$$C*
c$$$C*	  Ch(X,chi0) = - Q * exp(X) * d/dQ [ Z(Q) ]
c$$$C*
c$$$C*	with Q = X*sin(chi0), s0 = 1/sin(chi0).  Differentiating 
c$$$C*	inside the integral, we find that
c$$$C*
c$$$C*	  Z"(Q) + Z'(Q)/Q - Z(Q) = sqrt(s0**2-1) * exp(-Q*s0) / Q
c$$$C*
c$$$C*	giving us an inhomogeneous modified Bessel's equation of order 
c$$$C*	zero.  Following Rabenstein [Ra66, pp43-45,149] the solution 
c$$$C*	of this equation can be written as
c$$$C*
c$$$C*	  Z(Q) = A*I0(Q) + B*K0(Q) - sqrt(s0**2-1) 
c$$$C*	         * int{t=Q,infinity} [ exp(-t*s0) 
c$$$C*		   * ( I0(Q)*K0(t) - I0(t)*K0(Q) ) dt ] 
c$$$C*
c$$$C*	with coefficients A and B to be determined by matching 
c$$$C*	boundary conditions.
c$$$C*
c$$$C*	    Differentiating with respect to Q we obtain
c$$$C*
c$$$C*	  Ch(X,chi0) = X*sin(chi0)*exp(X)*( 
c$$$C*		- A*I1(X*sin(chi0)) + B*K1(X*sin(chi0)) 
c$$$C*		+ cos(chi0) * int{y=X,infinity} [ exp(-y) 
c$$$C*		  * ( I1(X*sin(chi0))*K0(y*sin(chi0))
c$$$C*		    + K1(X*sin(chi0))*I0(y*sin(chi0)) ) dy ] )
c$$$C*
c$$$C*	Applying the boundary condition Ch(X,0) = 1 requires that 
c$$$C*	B = 0.  Similarly, the requirement that Ch(X,chi0) approach 
c$$$C*	the finite value of sec(chi0) as X approaches infinity [Ch31b, 
c$$$C*	p486, Equation (12)] implies A = 0.  Thus we have
c$$$C*
c$$$C*	  Ch(X,chi0) = X*sin(chi0)*cos(chi0)*exp(X)*
c$$$C*		int{y=X,infinity} [ exp(-y) 
c$$$C*		  * ( I1(X*sin(chi0))*K0(y*sin(chi0))
c$$$C*		    + K1(X*sin(chi0))*I0(y*sin(chi0)) ) dy ]
c$$$C*
c$$$C*	The function atm8_chap_deq(X,chi0) evaluates this expression.
c$$$C*	Since explicit approximations are available for I1(z) and K1(z),
c$$$C*	the remaining challenge is evaluation of the integrals
c$$$C*
c$$$C*	  int{y=X,infinity} [ exp(-y) I0(y*sin(chi0)) dy ]
c$$$C*
c$$$C*	and
c$$$C*
c$$$C*	  int{y=X,infinity} [ exp(-y) K0(y*sin(chi0)) dy ]
c$$$C*
c$$$C*	which are accomplished by term-by-term integration of ascending
c$$$C*	and descending power series expansions of I0(z) and K0(z).
c$$$C*
c$$$C*  REFERENCES:
c$$$C*
c$$$C*	AS64	M. Abramowitz and I. A. Stegun, "Handbook of 
c$$$C*		Mathematical Functions," NBS AMS 55 (USGPO, 
c$$$C*		Washington, DC, June 1964, 9th printing, November 1970).
c$$$C*
c$$$C*	Ch31b	S. Chapman, "The Absorption and Dissociative or
c$$$C*		Ionizing Effect of Monochromatic Radiation in an
c$$$C*		Atmosphere on a Rotating Earth: Part II. Grazing
c$$$C*		Incidence," Proc. Phys. Soc. (London), _43_, 483-501 
c$$$C*		(1931).
c$$$C*
c$$$C*	Fi64	J. A. Fitzmaurice, "Simplfication of the Chapman
c$$$C*		Function for Atmospheric Attenuation," Appl. Opt. _3_,
c$$$C*		640 (1964).
c$$$C*
c$$$C*	Hi56a	F. B. Hildebrand, "Introduction to Numerical
c$$$C*		Analysis," (McGraw-Hill, New York, 1956).
c$$$C*
c$$$C*	Hu39	E. O. Hulburt, "The E Region of the Ionosphere," 
c$$$C*		Phys. Rev. _55_, 639-645 (1939).
c$$$C*
c$$$C*	PFT86	W. H. Press, B. P. Flannery, S. A. Teukolsky, and 
c$$$C*		W. T. Vetterling, "Numerical Recipes," (Cambridge, 
c$$$C*		1986).
c$$$C*
c$$$C*	Ra66	A. L. Rabenstein, "Introduction to Ordinary
c$$$C*		Differential Equations," (Academic, NY, 1966).
c$$$C*
c$$$C*	Re89	M. H. Rees, "Physics and Chemistry of the Upper
c$$$C*		Atmosphere," (Cambridge, 1989).
c$$$C*
c$$$C*	SG69	W. Swider, Jr., and M. E. Gardner, "On the Accuracy 
c$$$C*		of Chapman Function Approximations," Appl. Opt. _8_,
c$$$C*		725 (1969).
c$$$C*
c$$$C*	Sw64	W. Swider, Jr., "The Determination of the Optical 
c$$$C*		Depth at Large Solar Zenith Angles," Planet. Space 
c$$$C*		Sci. _12_, 761-782 (1964).
c$$$C
c$$$C  ####################################################################
c$$$C
c$$$C	Chapman function calculated by various methods
c$$$C
c$$$C	  Ch(X,chi0) = atm_chapman(X,chi0)   : real*4 entry
c$$$C	  Ch(X,chi0) = atm8_chapman(X,chi0)  : real*8 entry
c$$$C
c$$$C	Internal service routines - user should not call, except for
c$$$C	testing.
c$$$C
c$$$C	  Ch(X,chi0) = atm8_chap_asy(X,chi0) : asymptotic expansion
c$$$C	  Ch(X,chi0) = atm8_chap_deq(X,chi0) : differential equation
c$$$C	  Ch(X,chi0) = atm_chap_num(X,chi0)  : real*4 numerical integral
c$$$C	  Ch(X,chi0) = atm8_chap_num(X,chi0) : real*8 numerical integral
c$$$C
c$$$C  ####################################################################
c$$$
c$$$C  ====================================================================
c$$$C
c$$$C	These are the entries for the user to call.
c$$$C
c$$$C	chi0 can range from 0 to 180 in degrees.  For chi0 .gt. 90, the 
c$$$C	product X*(1-sin(chi0)) must not be too large, otherwise we 
c$$$C	will get an exponential overflow.
c$$$C
c$$$C	For chi0 .le. 90 degrees, X can range from 0 to thousands 
c$$$C	without overflow.
c$$$C
c$$$C  ====================================================================

      real*8 function atm_chapman( x, chi0 )
      real*8 atm8_chapman
      atm_chapman = atm8_chapman( dble(x), dble(chi0) )
      return
      end

c  ====================================================================

	real*8 function atm8_chapman( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)

	if( (x .le. 0) .or. (chi0 .le. 0) .or. (chi0 .ge. 180) ) then
	  atm8_chapman = 1
	  return
	end if

	if( chi0 .gt. 90 ) then
	  chi = 180 - chi0
	else
	  chi = chi0
	end if

	if( x .lt. 36 ) then
	  atm8_chapman = atm8_chap_deq(x,chi)
	else
	  atm8_chapman = atm8_chap_asy(x,chi)
	end if

	if( chi0 .gt. 90 ) then
	  atm8_chapman = 2*exp(x*2*sin((90-chi)/(2*rad))**2)
     .		* atm8_chap_xk1(x*sin(chi/rad)) - atm8_chapman
	end if

	return
	end

c$$$c  ====================================================================
c$$$c
c$$$c	this chapman function routine calculates
c$$$c
c$$$c	  ch(x,chi0) = atm8_chap_asy(x,chi0)
c$$$c		     = sum{n=0,3} [c(n) * int{t=t0,infinity} 
c$$$c			[ exp(-x*(t**2-t0**2) * (t**2-t0**2)**n dy ] ]
c$$$c
c$$$c	with t0**2 = 1 - sin(chi0)
c$$$c
c$$$c  ====================================================================

	real*8 function atm8_chap_asy( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	dimension c(0:3), xi(0:3), dn(0:3)
	common/atm8_chap_cm/fn(0:3)

	if( (x .le. 0) .or. (chi0 .le. 0) ) then
	  do i=0,3
	    fn(i) = 1
	  end do
	  go to 900
	end if

	sinchi = sin(chi0/rad)
	s1 = 1 + sinchi
	rx = sqrt(x)
	y0 = rx * sqrt( 2*sin( (90-chi0)/(2*rad) )**2 )

	c(0) = 1/sqrt(s1)
	fact = c(0)/s1
	c(1) = fact * (0.5d0+sinchi)
	fact = fact/s1
	c(2) = - fact * (0.125d0+0.5d0*sinchi)
	fact = fact/s1
	c(3) = fact * (0.0625d0+0.375d0*sinchi)

	call atm8_chap_gd3( y0, dn )
	fact = 2*rx
	do n=0,3
	  xi(n) = fact * dn(n)
	  fact = fact/x
	end do

	fn(0) = c(0) * xi(0)
	do i=1,3
	  fn(i) = fn(i-1) + c(i)*xi(i)
	end do

900	atm8_chap_asy = fn(3)
	return
	end

c  ====================================================================
c
c	this chapman function routine calculates
c
c	  ch(x,chi0) = atm8_chap_deq(x,chi0)
c		     = x * sin(chi0) * cos(chi0) * exp(x*sin(chi0))
c		       * int{y=x,infinity} [ exp(-y)*( 
c			 i1(x*sin(chi0))*k0(y*sin(chi0)) 
c			 + k1(x*sin(chi0))*i0(y*sin(chi0)) ) dy ]
c
c  ====================================================================

	real*8 function atm8_chap_deq( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	common/atm8_chap_cm/xi1,xk1,yi0,yk0

	if( (x .le. 0) .or. (chi0 .le. 0) ) go to 800
	alpha = x * sin(chi0/rad)

c  --------------------------------------------------------------------
c
c	this code fragment calculates
c
c	  yi0 = exp(x*(1-sin(chi0))) * cos(chi0) * 
c		int{y=x,infinity} [ exp(-y) * i0(y*sin(chi0)) dy ]
c
c  --------------------------------------------------------------------

	yi0 = atm8_chap_yi0( x, chi0 )

c  --------------------------------------------------------------------
c
c	this code fragment calculates
c
c	  yk0 = exp(x*(1+sin(chi0))) x * sin(chi0) * cos(chi0) * 
c		int{y=x,infinity} [ exp(-y) * k0(y*sin(chi0)) dy ]
c
c  --------------------------------------------------------------------

	yk0 = atm8_chap_yk0( x, chi0 )

c  --------------------------------------------------------------------
c
c	this code fragment calculates
c
c	  xi1 = exp(-x*sin(chi0)) * i1(x*sin(chi0))
c
c  --------------------------------------------------------------------

	xi1 = atm8_chap_xi1( alpha )

c  --------------------------------------------------------------------
c
c	this code fragment calculates
c
c	  xk1 = x*sin(chi0) * exp(x*sin(chi0)) * k1(x*sin(chi0))
c
c  --------------------------------------------------------------------

	xk1 = atm8_chap_xk1( alpha )

c  --------------------------------------------------------------------
c
c	combine the terms
c
c  --------------------------------------------------------------------

	atm8_chap_deq = xi1*yk0 + xk1*yi0
	go to 900

800	atm8_chap_deq = 1
900	return
	end

c  ====================================================================
c
c	this chapman function routine calculates
c
c	  ch(x,chi0) = atm_chap_num(x,chi0) = numerical integral
c
c  ====================================================================

	real*4 function atm_chap_num(x,chi0)
	implicit real*8(a-h,o-z)
	real*4 x, chi0
	parameter (rad=57.2957795130823208768d0)
	parameter (n=65,nfact=8)
	dimension factor(0:nfact)
	data factor/14175.0d0, 23552.0d0, -3712.0d0, 41984.0d0,
     .	  -18160.0d0, 41984.0d0, -3712.0d0, 23552.0d0, 7912.0d0/

	if( (chi0 .le. 0) .or. (chi0 .gt. 90) .or. (x .le. 0) ) then
	  atm_chap_num = 1
	  return
	end if

	x8 = x
	chi0rad = chi0/rad
	sinchi = sin(chi0rad)

	alpha0 = asin( (x8/(x8+28)) * sinchi )
	delta = (chi0rad - alpha0)/(n-1)

	sum = 0

	do i=1,n
	  alpha = -(i-1)*delta + chi0rad

	  if( (i .eq. 1) .or. (x .le. 0) ) then
	    f = 1/(1+cos(alpha))
	  else if( alpha .le. 0 ) then
	    f = 0
	  else
	    f = exp(-x8*(sinchi/sin(alpha)-1) ) /(1+cos(alpha))
	  end if

	  if( (i.eq.1) .or. (i.eq.n) ) then
	    fact = factor(nfact)/2
	  else
	    fact = factor( mod(i-2,nfact)+1 )
	  end if

	  sum = sum + fact*f
	end do

	atm_chap_num = 1 + x8*sinchi*sum*delta/factor(0)
	return
	end

c  ====================================================================
c
c	this chapman function routine calculates
c
c	  ch(x,chi0) = atm8_chap_num(x,chi0) = numerical integral
c
c  ====================================================================

	real*8 function atm8_chap_num(x,chi0)
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	parameter (n=601,nfact=8)
	dimension factor(0:nfact)
	data factor/14175.0d0, 23552.0d0, -3712.0d0, 41984.0d0,
     .	  -18160.0d0, 41984.0d0, -3712.0d0, 23552.0d0, 7912.0d0/

	if( (chi0 .le. 0) .or. (chi0 .gt. 90) .or. (x .le. 0) ) then
	  atm8_chap_num = 1
	  return
	end if

	chi0rad = chi0/rad
	sinchi = sin(chi0rad)

	alpha0 = asin( (x/(x+45)) * sinchi )
	delta = (chi0rad - alpha0)/(n-1)

	sum = 0

	do i=1,n
	  alpha = -(i-1)*delta + chi0rad

	  if( (i .eq. 1) .or. (x .le. 0) ) then
	    f = 1/(1+cos(alpha))
	  else if( alpha .le. 0 ) then
	    f = 0
	  else
	    f = exp(-x*(sinchi/sin(alpha)-1) ) /(1+cos(alpha))
	  end if

	  if( (i.eq.1) .or. (i.eq.n) ) then
	    fact = factor(nfact)/2
	  else
	    fact = factor( mod(i-2,nfact)+1 )
	  end if

	  sum = sum + fact*f
	end do

	atm8_chap_num = 1 + x*sinchi*sum*delta/factor(0)
	return
	end

c  ####################################################################
c
c	the following "bessel integral" routines return various 
c	combinations of integrals of bessel functions, powers, 
c	and exponentials, involving trigonometric functions of chi0.
c
c	for small values of z = x*sin(chi0) we expand
c
c	  i0(z) = sum{n=0,6} [ ai0(n) * z**(2*n) ]
c	  k0(z) = -log(z)*i0(z) + sum{n=0,6} [ ak0(n) * z**(2*n) ]
c
c	for large values of z we expand in reciprocal powers
c
c	  i0(z) = exp(z) * sum{n=0,8} [ bi0(n) * z**(-n-0.5) ]
c	  k0(z) = exp(-z) * sum{n=0,6} [ bk0(n) * z**(-n-0.5) ]
c
c	the expansion coefficients are calculated from those given 
c	by abramowitz and stegun [as64, pp378-9, section 9.8] and
c	press et al. [pft86, pp177-8, bessi0.for, bessk0.for].
c
c	for small values of x*sin(chi0) we break the integral
c	into two parts (with f(z) = i0(z) or k0(z)):
c
c	  int{y=x,infinity} [ exp(-y) * f(y*sin(chi0)) dy ]
c
c	    = int{y=x,x1} [ exp(-y) * f(y*sin(chi0)) dy ]
c	      + int{y=x1,infinity} [ exp(-y) * f(y*sin(chi0)) dy ]
c
c	where x1 = 3.75/sin(chi0) for i0 and 2/sin(chi0) for k0.
c
c	in the range y=x,x1 we integrate the term-by-term using
c
c	  int{z=a,b} [ exp(-z) * z**(2*n) dz ]
c	    = gamma(2*n+1,a) - gamma(2*n+1,b)
c
c	and a similar but more complicated formula for
c
c	  int{z=a,b} [ log(z) * exp(-z) * z**(2*n) dz ]
c
c	in the range y=x1,infinity we use
c
c	  int{z=b,infinity} [ exp(-z) * z**(-n-0.5) dz]
c	    = gamma(-n+0.5,b)
c
c  ####################################################################

c  ====================================================================
c
c	this bessel integral routine calculates
c
c	  yi0 = exp(x*(1-sin(chi0))) * cos(chi0) * 
c		int{y=x,infinity} [ exp(-y) * i0(y*sin(chi0)) dy ]
c
c  ====================================================================

	real*8 function atm8_chap_yi0( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	dimension qbeta(0:8), gg(0:6)
	dimension ai0(0:6), bi0(0:8)

        data ai0/ 1.0000000d+00, 2.4999985d-01, 1.5625190d-02,
     .      4.3393973d-04, 6.8012343d-06, 6.5601736d-08,
     .      5.9239791d-10/
        data bi0/ 3.9894228d-01, 4.9822200d-02, 3.1685484d-02,
     .     -8.3090918d-02, 1.8119815d+00,-1.5259477d+01,
     .      7.3292025d+01,-1.7182223d+02, 1.5344533d+02/

	theta = (90-chi0)/(2*rad)
	sint = sin(theta)
	cost = cos(theta)
	sinchi = sin(chi0/rad)
	coschi = cos(chi0/rad)
	sc1m = 2*sint**2	! = (1-sinchi)

	alpha = x * sinchi

	if( alpha .le. 0 ) then
	  atm8_chap_yi0 = 1
	else if( alpha .lt. 3.75d0 ) then
	  x1 = 3.75d0/sinchi
	  call atm8_chap_gg06( x, x1, gg )
	  if( x .le. 1 ) then
	    rho = 1
	  else
	    rho = 1/x
	  end if
	  f = (sinchi/rho)**2
	  sum = ai0(6)*gg(6)
	  do i=5,0,-1
	    sum = sum*f + ai0(i)*gg(i)
c	    write(*,1900)i,sum,gg(i)
c1900	format(i5,1p5d14.6)
	  end do
	  call atm8_chap_gq85( x1*sc1m, qbeta )
	  sum2 = bi0(8) * qbeta(8)
	  do n=7,0,-1
	    sum2 = sum2/3.75d0 + bi0(n)*qbeta(n)
	  end do
	  atm8_chap_yi0 = exp(-alpha)*coschi*sum 
     .		+ exp((x-x1)*sc1m)*sum2*cost*sqrt(2/sinchi)
	else
	  call atm8_chap_gq85( x*sc1m, qbeta )
	  sum = bi0(8) * qbeta(8)
	  do n=7,0,-1
	    sum = sum/alpha + bi0(n)*qbeta(n)
	  end do
	  atm8_chap_yi0 = sum * cost * sqrt( 2 / sinchi )
	end if
	return
	end

c  ====================================================================
c
c	this bessel integral routine calculates
c
c	  yk0 = exp(x*(1+sin(chi0))) x * sin(chi0) * cos(chi0) * 
c		int{y=x,infinity} [ exp(-y) * k0(y*sin(chi0)) dy ]
c
c  ====================================================================

	real*8 function atm8_chap_yk0( x, chi0 )
	implicit real*8(a-h,o-z)
	parameter (rad=57.2957795130823208768d0)
	dimension ai0(0:6), ak0(0:6), bk0(0:6)
	dimension gf(0:6), gg(0:6), qgamma(0:8)
	
        data ai0/ 1.0000000d+00, 2.4999985d-01, 1.5625190d-02,
     .      4.3393973d-04, 6.8012343d-06, 6.5601736d-08,
     .      5.9239791d-10/
        data ak0/ 1.1593152d-01, 2.7898274d-01, 2.5249154d-02,
     .      8.4587629d-04, 1.4975897d-05, 1.5045213d-07,
     .      2.2172596d-09/
        data bk0/ 1.2533141d+00,-1.5664716d-01, 8.7582720d-02,
     .     -8.4995680d-02, 9.4059520d-02,-8.0492800d-02,
     .      3.4053120d-02/

	theta = (90-chi0)/(2*rad)
	sint = sin(theta)
	cost = cos(theta)
	sinchi = sin(chi0/rad)
	sc1 = 1+sinchi
	coschi = sin(2*theta)

	alpha = x * sinchi
	gamma = x * sc1

	if( alpha .le. 0 ) then
	  atm8_chap_yk0 = 0
	else if( alpha .lt. 2 ) then
	  x1 = 2/sinchi
	  call atm8_chap_gfg06( x, x1, gf, gg )
	  if( x .le. 1 ) then
	    rho = 1
	  else
	    rho = 1/x
	  end if
	  sl = log(sinchi)
	  f = (sinchi/rho)**2
	  sum = -ai0(6)*gf(6) + (-sl*ai0(6)+ak0(6))*gg(6)
	  do i=5,0,-1
	    sum = sum*f - ai0(i)*gf(i) + (-sl*ai0(i)+ak0(i))*gg(i)
c	    write(*,1900)i,sum,gf(i),gg(i)
c1900	format(i5,1p5d14.6)
	  end do
	  call atm8_chap_gq85( x1*sc1, qgamma )
	  sum2 = bk0(6)*qgamma(6)
	  do i=5,0,-1
	    sum2 = sum2*0.5d0 + bk0(i)*qgamma(i)
c	    write(*,1900)i,sum2,bk0(i),qgamma(i)
	  end do
	  sum = sum + exp(x-x1-2)*sum2/sqrt(sinchi*sc1)
	  atm8_chap_yk0 = sum * exp(alpha) * alpha * coschi
	else
	  call atm8_chap_gq85( gamma, qgamma )
	  sum = bk0(6) * qgamma(6)
	  do i=5,0,-1
	    sum = sum/alpha + bk0(i)*qgamma(i)
	  end do
	  atm8_chap_yk0 = sum * sint * sqrt( 2 * sinchi ) * x
	end if

	return
	end

c  ####################################################################
c
c	the following "pure math" routines return various combinations
c	of bessel functions, powers, and exponentials.
c
c  ####################################################################

c  ====================================================================
c
c	this bessel function math routine returns
c
c	  xi1 = exp(-|z|) * i1(z)
c
c	following press et al [pft86, page 178, bessi1.for] and 
c	abrahamson and stegun [as64, page 378, 9.8.3, 9.8.4].
c
c  ====================================================================

	real*8 function atm8_chap_xi1( z )
	implicit real*8(a-h,o-z)
        dimension ai1(0:6), bi1(0:8)

        data ai1/ 5.00000000d-01, 6.2499978d-02, 2.6041897d-03,
     .      5.4244512d-05, 6.7986797d-07, 5.4830314d-09,
     .      4.1909957d-11/
        data bi1/ 3.98942280d-01,-1.4955090d-01,-5.0908781d-02,
     .      8.6379434d-02,-2.0399403d+00, 1.6929962d+01,
     .     -8.0516146d+01, 1.8642422d+02,-1.6427082d+02/

	if( z .lt. 0 ) then
	  az = -z
	else if( z .eq. 0 ) then
	  atm8_chap_xi1 = 0
	  return
	else
	  az = z
	end if
	if( az .lt. 3.75d0 ) then
	  z2 = z*z
	  sum = ai1(6)
	  do i=5,0,-1
	    sum = sum*z2 + ai1(i)
	  end do
	  atm8_chap_xi1 = z*exp(-az) * sum
	else
	  sum = bi1(8)
	  do i=7,0,-1
	    sum = sum/az + bi1(i)
	  end do
	  atm8_chap_xi1 = sum*sqrt(az)/z
	end if
	return
	end

c  ====================================================================
c
c	this bessel function math routine returns
c
c	  xk1 = z * exp(+z) * k1(z)
c
c	following press et al [pft86, page 179, bessk1.for] and 
c	abrahamson and stegun [as64, page 379, 9.8.7, 9.8.8].
c
c  ====================================================================

	real*8 function atm8_chap_xk1( z )
	implicit real*8(a-h,o-z)
        dimension ak1(0:6), bk1(0:6)

        data ak1/ 1.00000000d+00, 3.8607860d-02,-4.2049112d-02,
     .     -2.8370152d-03,-7.4976641d-05,-1.0781641d-06,
     .     -1.1440430d-08/
        data bk1/ 1.25331414d+00, 4.6997238d-01,-1.4622480d-01,
     .      1.2034144d-01,-1.2485648d-01, 1.0419648d-01,
     .     -4.3676800d-02/

	if( z .le. 0 ) then
	  atm8_chap_xk1 = 1
	else if( z .lt. 2 ) then
	  xz = exp(z)
	  z2 = z*z
	  sum = ak1(6)
	  do i=5,0,-1
	    sum = sum*z2 + ak1(i)
	  end do
	  atm8_chap_xk1 = xz * ( sum 
     *		+ z*log(z/2)*atm8_chap_xi1(z)*xz )
	else
	  sum = bk1(6)
	  do i=5,0,-1
	    sum = sum/z + bk1(i)
	  end do
	  atm8_chap_xk1 = sum*sqrt(z)
	end if

	return
	end

c  ####################################################################
c
c	the following "pure math" routines return various combinations
c	of the error function, powers, and exponentials.
c
c  ####################################################################

c  ====================================================================
c
c	this error function math routine returns
c
c	  xerfc(x) = exp(x**2)*erfc(x)
c
c	following press et al. [pft86, p164, erfcc.for]
c
c  ====================================================================

	real*8 function atm8_chap_xerfc(x)
	implicit real*8(a-h,o-z)
        t=1.0d0/(1.0d0+0.5d0*x)
	atm8_chap_xerfc =
     .	  t*exp( -1.26551223d0 +t*(1.00002368d0 +t*( .37409196d0
     .       +t*(  .09678418d0 +t*(-.18628806d0 +t*( .27886807d0
     .	     +t*(-1.13520398d0 +t*(1.48851587d0 +t*(-.82215223d0
     .       +t*   .17087277d0) ))))))))
        return
        end

c  ####################################################################
c
c	the following "pure math" routines return various combinations
c	of exponential integrals, powers, and exponentials.
c
c  ####################################################################

c  ====================================================================
c
c	this exponential math routine evaluates
c
c	  zxe1(x) = x*exp(x) int{y=1,infinity} [ exp(-x*y)/y dy ]
c
c	following abramowitz and stegun [as64, p229;231, equations
c	5.1.11 and 5.1.56]
c
c  ====================================================================

	real*8 function atm8_chap_zxe1(x)
	implicit real*8(a-h,o-z)
	parameter (gamma = 0.5772156649015328606d0)
	dimension ae1(0:4), be1(0:4), cein(1:10)

	data ae1/1.0d0, 8.5733287401d0, 18.0590169730d0,
     .	    8.6347608925d0, 0.2677737343d0 /
	data be1/1.0d0, 9.5733223454d0, 25.6329561486d0,
     .	    21.0996530827d0, 3.9584969228d0/
        data cein/ 1.00000000d+00,-2.50000000d-01, 5.55555556d-02,
     .    -1.0416666667d-02, 1.6666666667d-03,-2.3148148148d-04,
     .     2.8344671202d-05,-3.1001984127d-06, 3.0619243582d-07,
     .    -2.7557319224d-08/

	if( x .le. 0 ) then
	  atm8_chap_zxe1 = 0
	else if( x .le. 1 ) then
	  sum = cein(10)
	  do i=9,1,-1
	    sum = sum*x + cein(i)
	  end do
	  atm8_chap_zxe1 = x*exp(x)*( x * sum - log(x) - gamma )
	else
	  top = ae1(4)
	  bot = be1(4)
	  do i=3,0,-1
	    top = top/x + ae1(i)
	    bot = bot/x + be1(i)
	  end do
	  atm8_chap_zxe1 = top/bot
	end if
	return
	end

c  ####################################################################
c
c	the following "pure math" routines return various combinations
c	of incomplete gamma functions, powers, and exponentials.
c
c  ####################################################################

c  ====================================================================
c
c	this gamma function math routine calculates
c
c	dn(n) = int{t=z,infinity}
c		[ exp( -(t**2-z**2) ) * (t**2-z**2)**n dt ]
c
c  ====================================================================

	subroutine atm8_chap_gd3( z, dn )
	implicit real*8(a-h,o-z)
	parameter (rpi=1.7724538509055160273d0)
	dimension dn(0:3), xg(0:3)

	if( z .le. 0 ) then
	  dn(0) = rpi/2
	  do i=1,3
	    dn(i) = (i-0.5d0)*dn(i-1)
	  end do
	  return
	end if

	z2 = z*z
	if( z .ge. 7 ) r = 1/z2

	if( z .lt. 14 ) then
	  z4 = z2*z2
	  xg(0) = rpi * atm8_chap_xerfc(z)
	  xg(1) = 0.5d0*xg(0) + z
	  xg(2) = 1.5d0*xg(1) + z*z2
	  dn(0) = 0.5d0*xg(0)
	  dn(1) = 0.5d0*(xg(1)-z2*xg(0))
	  dn(2) = 0.5d0*(xg(2)-2*z2*xg(1)+z4*xg(0))
	else
	  dn(0) = ( 1 + r*(-0.5d0 +r*(0.75d0 +r*(-1.875d0
     .		+r*6.5625d0) ) ) )/(2*z)
	  dn(1) = ( 1 + r*(-1.0d0 +r*(2.25d0 +r*(-7.5d0
     .		+r*32.8125d0) ) ) )/(2*z)
	  dn(2) = ( 2 + r*(-3.0d0 +r*(9.00d0 +r*(-37.5d0
     .		+r*196.875d0) ) ) )/(2*z)
	end if

	if( z .lt. 7 ) then
	  z6 = z4*z2
	  xg(3) = 2.5d0*xg(2) + z*z4
	  dn(3) = 0.5d0*(xg(3)-3*z2*xg(2)+3*z4*xg(1)-z6*xg(0))
	else
	  dn(3) = ( 6 + r*(-12.0d0 +r*(45.0d0 +r*(-225.0d0
     .		+r*1378.125d0) ) ) )/(2*z)
	end if

	return
	end

c  ====================================================================
c
c	this gamma function math routine calculates
c
c	  gf06(n) = g(n,x) * int{y=x,z} [log(y) * exp(-y) * y**(2*n) dy]
c
c	and
c
c	  gg06(n) = g(n,x) * int{y=x,z} [ exp(-y) * y**(2*n) dy ]
c	          = g(n,x) * ( gamma(2*n+1,x) - gamma(2*n+1,z) )
c
c	for n=0,6, with g(n,x) = exp(x) * max(1,x)**(-2*n)
c
c  ====================================================================

	subroutine atm8_chap_gfg06( x, z, gf06, gg06 )
	implicit real*8 (a-h,o-z)
	parameter (gamma = 0.5772156649015328606d0)
	dimension gf06(0:6), gg06(0:6)
	dimension gh13x(13), gh13z(13), rgn(13), delta(13)
	call atm8_chap_gh13( x, x, gh13x )
	call atm8_chap_gh13( x, z, gh13z )
	if( x .le. 1 ) then
	  rho = 1
	else
	  rho = 1/x
	end if

	delta(1) = 0
	delta(2) = ( gh13x(1) - gh13z(1) ) * rho
	rgn(1) = 1
	rgn(2) = rho
	do n=2,12
	  delta(n+1) = rho*( n*delta(n) + gh13x(n) - gh13z(n) )
	  rgn(n+1) = (n*rho)*rgn(n)

	end do

	if( x .gt. 0 ) then
	  xe1_x = atm8_chap_zxe1(x)/x
	  xlog = log(x)
	end if
	if( z .gt. 0 ) then
	  xe1_z = exp(x-z)*atm8_chap_zxe1(z)/z
	  zlog = log(z)
	end if

	do k=0,6
	  n = 2*k+1
	  if( x .le. 0 ) then
	    gf06(k) = -gamma*rgn(n) + delta(n)
	  else
	    gf06(k) = xlog*gh13x(n) + rgn(n)*xe1_x + delta(n)
	  end if
	  if( z .le. 0 ) then
	    gf06(k) = gf06(k) + gamma*rgn(n)
	  else
	    gf06(k) = gf06(k) - (zlog*gh13z(n) + rgn(n)*xe1_z)
	  end if
	  gg06(k) = gh13x(n) - gh13z(n)
	end do

	return
	end

c  ====================================================================
c
c	this gamma function math routine calculates
c
c	  gg06(n) = g(n,x) * int{y=x,z} [ exp(-y) * y**(2*n) dy ]
c	          = g(n,x) * ( gamma(2*n+1,x) - gamma(2*n+1,z) )
c
c	for n=0,6, with g(n,x) = exp(x) * max(1,x)**(-2*n)
c
c  ====================================================================

	subroutine atm8_chap_gg06( x, z, gg06 )
	implicit real*8 (a-h,o-z)
	dimension gg06(0:6), gh13x(13), gh13z(13)
	call atm8_chap_gh13( x, x, gh13x )
	call atm8_chap_gh13( x, z, gh13z )
	do n=0,6
	  gg06(n) = gh13x(2*n+1) - gh13z(2*n+1)
	end do
	return
	end

c  ====================================================================
c
c	this gamma function math routine calculates
c
c	  gh13(n) = f(n,x) * int{y=z,infinity} [exp(-y) * y**(n-1) dy]
c	          = f(n,x) * gamma(n,z)
c
c	for n=1,13, with f(n,x) = exp(x) * max(1,x)**(-n+1)
c
c  ====================================================================

	subroutine atm8_chap_gh13( x, z, gh13 )
	implicit real*8 (a-h,o-z)
	dimension gh13(13), tab(12)

	if( z .le. 0 ) then
	  gh13(1) = 1
	  do n=1,12
	    gh13(n+1) = n*gh13(n)
	  end do
	  return
	end if

	if( x .le. 1 ) then
	  rho = 1
	else
	  rho = 1/x
	end if
	rhoz = rho * z
	exz = exp(x-z)
	tab(12) = exp( (x-z) + 12*log(rhoz) )
	do n=11,1,-1
	  tab(n) = tab(n+1)/rhoz
	end do
	gh13(1) = exz
	do n=1,12
	  gh13(n+1) = rho*n*gh13(n) + tab(n)
	end do
	return
	end

c  ====================================================================
c
c	this gamma function math subroutine calculates
c
c	  qn(x) = x**n * exp(x) * gamma(-n+0.5,x), n=0,8
c	    = x**n * exp(x) * int{y=x,infinity} [exp(-y)*y**(-n-0.5)dy]
c
c	for x .lt. 2 we first calculate
c
c	  q0(x) = sqrt(pi)*exp(x)*erfc(sqrt(x)) = exp(x)*gamma(0.5,x)
c
c	and use upward recursion.  else, we first calculate
c
c	  q8(x) = x**8 * exp(x) * gamma(-7.5,x)
c
c	following press et al. [pft86, pp162-63, gcf.for] and then
c	recur downward.  also see abramowitz and stegun [as64, 6.5].
c
c  ====================================================================

	subroutine atm8_chap_gq85( x, qn )
	implicit real*8(a-h,o-z)
	parameter (rpi=1.7724538509055160273d0)
	parameter (itmax=100,eps=3.0d-9)
	dimension qn(0:8)

	if( x .le. 0 ) then
	  qn(0) = rpi
	  do i=1,8
	    qn(i) = 0
	  end do
	  return
	end if

	rx = sqrt(x)

	if( x .lt. 2 ) then
	  qn(0) = rpi * atm8_chap_xerfc( rx )
	  do n=1,8
	    qn(n) = ( -rx*qn(n-1) + 1 ) * rx / ( n - 0.5d0 )
	  end do
	else
          gold=0.0d0
	  a0=1.0d0
	  a1=x
	  b0=0.0d0
	  b1=1.0d0
	  fac=1.0d0
	  do 11 n=1,itmax
	    an= (n)
	    ana=an + 7.5d0
	    a0=(a1+a0*ana)*fac
	    b0=(b1+b0*ana)*fac
	    anf=an*fac
	    a1=x*a0+anf*a1
	    b1=x*b0+anf*b1
	    fac=1./a1
            g=b1*fac
	    test = g*eps
	    del = g - gold
	    if( test .lt. 0 ) test = - test
	    if( (del .ge. -test) .and. (del .le. test) ) go to 12
	    gold=g
11        continue
12	  qn(8) = g * rx
	  do n=8,1,-1
	    qn(n-1) = ( (-n+0.5d0)*qn(n)/rx + 1 ) / rx
	  end do
	end if

	return
	end


