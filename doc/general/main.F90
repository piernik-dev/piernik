!> \mainpage PIERNIK DOCUMENTATION (UNDER CONSTRUCTION)
!!
!!\par Introduction
!!
!!PIERNIK is an MHD code created at Centre for Astronomy, Nicolaus Copernicus University in Torun, Poland. Current version of the code uses a simple, conservative numerical scheme, which is known as Relaxing TVD scheme (RTVD). General mathematical context of the relaxation and relaxing systems of hyperbolic conservation laws, and related numerical schemes, was presented by Jin & Xin (1995). A particular realization of the Relaxing TVD was developed by Trac & Pen (2003) and Pen et al. (2003), who presented the numerical method in a pedagogical way, and provided short, publicly available HD and MHD codes. These codes rely on a dimensionally split, second order algorithm in  space and time.  The Relaxing TVD scheme is easily extendible to account for additional fluid components: multiple fluids, dust, cosmic rays, and additional physical processes, such as fluid interactions and Ohmic resistivity  effects.  The simplicity and a small number of floating point operations of the basic  algorithm is reflected in a performance of $10^5$ zone-cycles/s (on  single-core 2 GHz processors). 
!!
!!
!! \par PIERNIK source code and documentation 
!!
!!Public version of PIERNIK code is available via the web-page: \url{http://piernik.astri.uni.torun.pl}. The web-page informs how to access the source code,  which is maintained in an SVN repository, together with on-line documentation describing details of code utilization. An overview of PIERNIK capabilities is in conference papers (Hanasz et al. 2008a,b,c and 2009a, see: \url{http://arxiv.org/abs/0812.2161}, \url{http://arxiv.org/abs/0812.2799}, \url{http://arxiv.org/abs/0812.4839}, \url{http://arxiv.org/abs/0901.0104}), and is also available in doc/papers subdirectory. Lectures notes on PIERNIK (slides in pdf format, as well as two sections prepared for Lecture Notes in Physics, Springer (in press)) are placed in doc/evora subdirectory. 
!!
!<

