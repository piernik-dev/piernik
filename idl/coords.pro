PRO COORDS, coord_sys

COMMON coordsys    , coordsys, coordnames $
                   , coord1,coord2,coord3,plane1,plane2,plane3

coordsys = coord_sys

xyz = ['x','y','z']
zrp = ['z','r','p']
rtp = ['r','t','p']

IF(coordsys EQ 'xyz') THEN BEGIN
  coordnames=xyz
  coord1 = 'x'
  plane1 = 'yz'
  coord2 = 'y'
  plane2 = 'xz'
  coord3 = 'z'
  plane3 = 'xy'

ENDIF ELSE IF(coordsys EQ 'zrp') THEN BEGIN
  coordnames=zrp
  coord1 = 'z'
  plane1 = 'rp'
  coord2 = 'r'
  plane2 = 'pz'
  coord3 = 'p'
  plane3 = 'zr'
ENDIF ELSE IF(coordsys EQ 'rtp') THEN BEGIN
  coordnames=rtp
  coord1 = 'r'
  plane1 = 'tp'
  coord2 = 't'
  plane2 = 'pr'
  coord3 = 'p'
  plane3 = 'rt'
ENDIF ELSE  BEGIN
  coordnames=[' ',' ',' ']
ENDELSE

END
