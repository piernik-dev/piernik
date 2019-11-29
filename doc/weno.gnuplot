#!/usr/bin/gnuplot -persist

# Note_0: Equation numbers according to Yamaleev, Carpenter, Journal of Computational Physics 228 (2009) 3025–3047

# Let's focus on i-th cell. We can shift our data vector in such a way:
#   u_i = 0
#   u_{i-1} = x = r * cos(phi)
#   u_(i+1) = y = r * sin(phi)
# Now, for simplicity we may put r = 1

# We're looking for the value at a point placed in a middle point between u_i and
# u_{i+1}, let's call it f = u_{i+1/2}

# Simple Linear interpolation to a middle point between u_i and u_{i+1}
# gives f = y/2 = sin(phi)/2

# An extrapolation to the same point from u_i and u_{i-1} gives
# f = -x/2 = -cos(phi)/2

# Any sane interpolation based on the points u_{i-1}, u_{i} and u_{i+1}
# should return values between -x/2, y/2 (inclusive).

# We can also fit a parabola to these three points and get
# f = (3y - x)/8 = (3sin(phi) - cos(phi)/8

# Now let's calculate beta_0 and beta_1 (Eq. 20) and tau (Eq. 22)
b0(x,y) = y**2
b1(x,y) = x**2
tau(x,y) = (x+y)**2

# We have two prescriptions for weights.
# Classic WENO (alpha_0 and alpha_1, Eq. 19)
# with e = epsilon (numerical safety factor)
wa0(x,y,e) = 2. / 3. / (e + b0(x,y))**2
wa1(x,y,e) = 1. / 3. / (e + b1(x,y))**2

# Or ESWENO (alpha_0 and alpha_1, Eq. 21)
sa0(x,y,e) = 2. / 3. * (1 + tau(x,y) / (e + b0(x,y)))
sa1(x,y,e) = 1. / 3. * (1 + tau(x,y) / (e + b1(x,y)))

# Now we can calculate w^0_{j+1/2} (Eq. 18) for both cases (WENO and ESWENO,
# respectively). Note that w^1 = 1 - w^0 and we don't have to define it explicitly:
ww0(x,y,e) = wa0(x,y,e) / (wa0(x,y,e) + wa1(x,y,e))
sw0(x,y,e) = sa0(x,y,e) / (sa0(x,y,e) + sa1(x,y,e))

# Finally we can evaluate y or f^W_{j+1/2} (Eq. 14 and 15, WENO and ESWENO,
# respectively)
wf(x,y,e)=ww0(x,y,e)*y/2 - (1-ww0(x,y,e))*x/2
sf(x,y,e)=sw0(x,y,e)*y/2 - (1-sw0(x,y,e))*x/2

# Now we can plot:
# * linear interpolation
# * linear extrapolation
# * parabolic interpolation
# * ESWENO with vanishing epsilon
# * WENO  with vanishing epsilon
# * (ES)WENO dominated by epsilon (e.g. data vector consists of tiny values)

set xlabel "atan2(y,x)/pi"
set key bottom
set ytics .1
set xtics .25
set grid
set samples 1000
set title "See the comments in the file for detailed explanations."
plot [-1:1] sin(pi*x)/2 title "linear interpolation", \
  cos(pi*x)/(-2) title "linear extrapolation", \
  (3*sin(pi*x)-cos(pi*x))/8 title "parabolic interpolation", \
  sf(cos(pi*x), sin(pi*x), 0.) title "ESWENO epsilon = 0", \
  wf(cos(pi*x), sin(pi*x), 0.) title "WENO epsilon = 0", \
  sf(cos(pi*x), sin(pi*x), 1000.) title "(ES)WENO epsilon = 1000"

# Interpretation:
# * Points at angles -pi/4 and 3pi/4 represent reconstruction of linear
#   slope - all schemes should cross them. The slope there comes from the
#   choice of weights (2/3 and 1/3 for (ES)WENO, Eq. 20)
# * (ES)WENO returning 0 for angles 0, +/-pi/2 and +/-pi represent the
#   property that tries to preserve sharp steps by extending the
#   interpolation from the less sloppy pair of points - that's why (ES)WENO
#   curves touch linear "interpolation" and "extrapolation" curves.
# * The value at angle pi/4 corresponds to reconstruction of local minimum
#   (u_j < min(u_{j-1}, u_{j+1}) and at -3pi/4 we have a local maximum
#   (u_j > max(u_{j-1}, u_{j+1}).
# * When epsilon is chosen too big, the scheme degrades to something close
#   to parabolic interpolation and cannot recontruct sharp steps.
