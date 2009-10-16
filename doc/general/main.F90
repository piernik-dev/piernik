!> \mainpage PIERNIK 
!>
!>\section Introduction
!!
!!PIERNIK is an MHD code (see Hanasz et al. \cite{hanasz-etal-09a,hanasz-etal-09b,hanasz-etal-09c,hanasz-etal-09d})
!!created at Centre for Astronomy, Nicolaus Copernicus University in Toru\'n, Poland. Current version of the code 
!!uses a simple, conservative numerical scheme, which is known as Relaxing TVD scheme (RTVD). General mathematical 
!!context of the relaxation and relaxing systems of hyperbolic conservation laws, and related numerical schemes, 
!!was presented by Jin \& Xin~\cite{jin-xin-95}. A particular realization of the Relaxing TVD was developed by 
!!Trac~\&~Pen~\cite{trac-pen-03} and Pen \etal~\cite{pen-etal-03}, who presented the numerical method in a pedagogical 
!!way, and provided short, publicly available HD and MHD codes. These codes rely on a dimensionally split, second order 
!!algorithm in  space and time.  The Relaxing TVD scheme is easily extendible to account for additional fluid components:
!! multiple fluids, dust, cosmic rays, and additional physical processes, such as fluid interactions and ohmic resistivity 
!!effects.  The simplicity and a small number of floating point operations of the basic algorithm is reflected in a 
!!performance of $10^5$ zone--cycles/s (on  single--core 2 GHz processors).
!!
!!Public version of PIERNIK code is available via the web-page:
!!\url{http://piernik.astri.uni.torun.pl}. The web--page informs how to access the source code,  which  is
!!maintained in an SVN repository, together with on--line documentation describing details of code utilization.

