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
!!
!!\section The Relaxing TVD scheme for MHD simulations
!!
!!Numerical algorithm of PIERNIK code is based on  the conservative form of the MHD system equations
!!\f[
!!  \partial _t \vec{u} + \partial_x \vec{F}{(\vec{u},\vec{B})} + \partial_y \vec{G}{(\vec{u},\vec{B})} + \partial_z \vec{H}{(\vec{u},\vec{B})} = \vec{S}(\vec{u}).
!!\f]
!!where \f$\vec{u} = \left( \rho, m_x, m_y, m_z, e \right)^T\f$
!!is a vector of conservative fluid variables: gas density, three components of momentum and total energy density, respectively, while 
!!\f$\vec{F}{(\vec{u},\vec{B})}\f$, \f$\vec{G}{(\vec{u},\vec{B})}\f$, \f$\vec{H}{(\vec{u},\vec{B})}\f$ are fluxes of the fluid variables in x, y and z directions, respectively, and \f$\vec{S}(\vec{u})\f$ is a function representing source terms.
!!The flux functions are given by 
!!\f[
!!  \vec{F}{(\vec{u})} = 
!!  \left(\begin{array}{c}
!!    \rho v_x \\
!!    \rho v_x^2 + p_* - B_x^2 \\
!!    \rho v_x v_y - B_x B_y\\
!!    \rho v_x v_z - B_x B_z\\
!!    (e + p)v_x - \vec{B} \cdot \vec{v} \; B_x 
!!  \end{array}\right),
!!\f]
!!\f[
!!  \vec{G}{(\vec{u})} = 
!!  \left(\begin{array}{c}
!!    \rho v_y \\
!!    \rho v_y v_x  - B_y B_x\\
!!    \rho v_y^2 + p_* - B_y^2 \\
!!    \rho v_y v_z  - B_y B_z\\
!!    (e + p)v_y - \vec{B} \cdot \vec{v} \; B_y 
!!  \end{array}\right),
!!\f]
!!\f[
!!  \vec{H}{(\vec{u})} = 
!!  \left(\begin{array}{c}
!!    \rho v_z \\
!!    \rho v_z v_x  - B_z B_x\\
!!    \rho v_z v_y  - B_z B_y\\
!!    \rho v_z^2 + p_* - B_z^2 \\
!!    (e + p)v_z - \vec{B} \cdot \vec{v} \; B_z 
!!  \end{array}\right),
!!\f]
!!where \f$p_* = p + B^2/2\f$, \f$e= e_{th} + \frac{1}{2} \rho v^2 + B^2/2\f$,  are the total pressure and total energy density, 
!!while \f$e_{th}\f$ is thermal energy density and  \f$e_{mag} = B^2/2\f$ is the magnetic energy density.
!!\image html tvd1l.png "Caption"
!!\image latex tvd1l.eps "My application"
!!\image html tvd1r.png
!!\image latex tvd1r.eps "My application"

