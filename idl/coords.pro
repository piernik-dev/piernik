PRO COORDS, crdsys, cs

xyz = ['x','y','z']
zrp = ['z','r','p']
rtp = ['r','t','p']

IF(crdsys EQ 'xyz') THEN BEGIN
  crdnam =  xyz
  coord1 = 'x'
  plane1 = 'yz'
  coord2 = 'y'
  plane2 = 'xz'
  coord3 = 'z'
  plane3 = 'xy'

ENDIF ELSE IF(crdsys EQ 'zrp') THEN BEGIN
  crdnam =  zrp
  coord1 = 'z'
  plane1 = 'rp'
  coord2 = 'r'
  plane2 = 'pz'
  coord3 = 'p'
  plane3 = 'zr'
ENDIF ELSE IF(crdsys EQ 'rtp') THEN BEGIN
  crdnam =  rtp
  coord1 = 'r'
  plane1 = 'tp'
  coord2 = 't'
  plane2 = 'pr'
  coord3 = 'p'
  plane3 = 'rt'
ENDIF ELSE  BEGIN
  coordnames=[' ',' ',' ']
ENDELSE

cs={crdsys:crdsys, $
    crdnam:crdnam, $
    coord1:coord1, $
    coord2:coord2, $
    coord3:coord3, $
    plane1:plane1, $
    plane2:plane2, $
    plane3:plane3}
END
