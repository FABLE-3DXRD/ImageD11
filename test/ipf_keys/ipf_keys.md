# MTEX 7.2.0 direction colour keys, lattice terms

Normals (hkl), directions [uvw]; ~ marks a direction with irrational indices for this cell.
A view is out-of-screen|east. Only views whose sector is a fundamental region are listed.
hue0 is the direction hue 0 points to from the centre. Colours compared on 500 random directions.

## 1

Laue -1, 1 operations, test cell [5, 6, 7] [80, 95, 105]

1. **TSL** N (0 0 1) | vertices - | center ~[0.0633 -0.1889 1.0000] | ref ~[1.0000 0.2094 0.0311] | hue0 ref  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 6.1e-14)
2. **TSL** N (1 0 0) | vertices - | center ~[1.0000 0.2094 0.0311] | ref ~[1.0000 0.2094 0.0311] | hue0 ref  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **TSL** N (-1 0 0) | vertices - | center ~[-1.0000 -0.2094 -0.0311] | ref ~[1.0000 0.2094 0.0311] | hue0 ref  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **TSL** N (0 1 0) | vertices - | center ~[0.2947 1.0000 -0.1305] | ref ~[1.0000 0.2094 0.0311] | hue0 ref  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
5. **TSL** N (0 -1 0) | vertices - | center ~[-0.2947 -1.0000 0.1305] | ref ~[1.0000 0.2094 0.0311] | hue0 ref  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
6. **TSL** N (0 0 -1) | vertices - | center ~[-0.0633 0.1889 -1.0000] | ref ~[1.0000 0.2094 0.0311] | hue0 ref  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix C1 TSL** N ~(-0.0623 0.1488 1.0000) | vertices - | center [0 0 1]  
    nearest MTEX key 1, max |dRGB| 0.995 on 497/500 directions > 0.01
1. **HSV** N - | refl (1 0 0) | vertices - | whiteCenter ~[1.0000 0.2094 0.0311] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b
2. **HSV** N - | refl (1 0 0) | vertices - | whiteCenter ~[1.0000 0.2094 0.0311] | ref [0 1 0] | hue0 ref  
   views: +a*|+b -a*|+b
3. **HSV** N - | refl (1 0 0) | vertices - | whiteCenter ~[1.0000 0.2094 0.0311] | ref [0 -1 0] | hue0 ref  
   views: +a*|-b -a*|-b
4. **HSV** N - | refl (1 0 0) | vertices - | whiteCenter ~[1.0000 0.2094 0.0311] | ref [0 0 1] | hue0 ref  
   views: +a*|+c -a*|+c
5. **HSV** N - | refl (1 0 0) | vertices - | whiteCenter ~[1.0000 0.2094 0.0311] | ref [0 0 -1] | hue0 ref  
   views: +a*|-c -a*|-c
6. **HSV** N - | refl (1 0 0) | vertices - | whiteCenter ~[1.0000 0.2094 0.0311] | ref [1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
7. **HSV** N - | refl (1 0 0) | vertices - | whiteCenter ~[1.0000 0.2094 0.0311] | ref [1 0 0] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
8. **HSV** N - | refl (1 0 0) | vertices - | whiteCenter ~[1.0000 0.2094 0.0311] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **pymtex HSV** N - | vertices - | center [0 0 1]  
    nearest MTEX key 6, max |dRGB| 0.888 on 500/500 directions > 0.01

## -1

Laue -1, 2 operations, test cell [5, 6, 7] [80, 95, 105]

1. **TSL** N (0 0 1) | vertices - | center ~[0.0633 -0.1889 1.0000] | ref ~[1.0000 0.2094 0.0311] | hue0 ref  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 6.1e-14)
2. **TSL** N (1 0 0) | vertices - | center ~[1.0000 0.2094 0.0311] | ref ~[1.0000 0.2094 0.0311] | hue0 ref  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **TSL** N (-1 0 0) | vertices - | center ~[-1.0000 -0.2094 -0.0311] | ref ~[1.0000 0.2094 0.0311] | hue0 ref  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **TSL** N (0 1 0) | vertices - | center ~[0.2947 1.0000 -0.1305] | ref ~[1.0000 0.2094 0.0311] | hue0 ref  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
5. **TSL** N (0 -1 0) | vertices - | center ~[-0.2947 -1.0000 0.1305] | ref ~[1.0000 0.2094 0.0311] | hue0 ref  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
6. **TSL** N (0 0 -1) | vertices - | center ~[-0.0633 0.1889 -1.0000] | ref ~[1.0000 0.2094 0.0311] | hue0 ref  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix Ci TSL** N ~(-0.0623 0.1488 1.0000) | vertices - | center [0 0 1]  
    nearest MTEX key 1, max |dRGB| 0.995 on 497/500 directions > 0.01
1. **HSV** N (0 0 1) | refl - | vertices - | whiteCenter ~[0.0633 -0.1889 1.0000] | ref [0 1 0] | hue0 ref  
   views: default +c*|+b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (1 0 0) | refl - | vertices - | whiteCenter ~[1.0000 0.2094 0.0311] | ref [0 1 0] | hue0 ref  
   views: +a*|+b
3. **HSV** N (1 0 0) | refl - | vertices - | whiteCenter ~[1.0000 0.2094 0.0311] | ref [0 -1 0] | hue0 ref  
   views: +a*|-b
4. **HSV** N (-1 0 0) | refl - | vertices - | whiteCenter ~[-1.0000 -0.2094 -0.0311] | ref [0 1 0] | hue0 ref  
   views: -a*|+b
5. **HSV** N (-1 0 0) | refl - | vertices - | whiteCenter ~[-1.0000 -0.2094 -0.0311] | ref [0 -1 0] | hue0 ref  
   views: -a*|-b
6. **HSV** N (1 0 0) | refl - | vertices - | whiteCenter ~[1.0000 0.2094 0.0311] | ref [0 0 1] | hue0 ref  
   views: +a*|+c
7. **HSV** N (1 0 0) | refl - | vertices - | whiteCenter ~[1.0000 0.2094 0.0311] | ref [0 0 -1] | hue0 ref  
   views: +a*|-c
8. **HSV** N (-1 0 0) | refl - | vertices - | whiteCenter ~[-1.0000 -0.2094 -0.0311] | ref [0 0 1] | hue0 ref  
   views: -a*|+c
9. **HSV** N (-1 0 0) | refl - | vertices - | whiteCenter ~[-1.0000 -0.2094 -0.0311] | ref [0 0 -1] | hue0 ref  
   views: -a*|-c
10. **HSV** N (0 1 0) | refl - | vertices - | whiteCenter ~[0.2947 1.0000 -0.1305] | ref [1 0 0] | hue0 ref  
   views: +b*|+a
11. **HSV** N (0 1 0) | refl - | vertices - | whiteCenter ~[0.2947 1.0000 -0.1305] | ref [-1 0 0] | hue0 ref  
   views: +b*|-a
12. **HSV** N (0 -1 0) | refl - | vertices - | whiteCenter ~[-0.2947 -1.0000 0.1305] | ref [1 0 0] | hue0 ref  
   views: -b*|+a
13. **HSV** N (0 -1 0) | refl - | vertices - | whiteCenter ~[-0.2947 -1.0000 0.1305] | ref [-1 0 0] | hue0 ref  
   views: -b*|-a
14. **HSV** N (0 1 0) | refl - | vertices - | whiteCenter ~[0.2947 1.0000 -0.1305] | ref [0 0 1] | hue0 ref  
   views: +b*|+c
15. **HSV** N (0 1 0) | refl - | vertices - | whiteCenter ~[0.2947 1.0000 -0.1305] | ref [0 0 -1] | hue0 ref  
   views: +b*|-c
16. **HSV** N (0 -1 0) | refl - | vertices - | whiteCenter ~[-0.2947 -1.0000 0.1305] | ref [0 0 1] | hue0 ref  
   views: -b*|+c
17. **HSV** N (0 -1 0) | refl - | vertices - | whiteCenter ~[-0.2947 -1.0000 0.1305] | ref [0 0 -1] | hue0 ref  
   views: -b*|-c
18. **HSV** N (0 0 1) | refl - | vertices - | whiteCenter ~[0.0633 -0.1889 1.0000] | ref [1 0 0] | hue0 ref  
   views: +c*|+a
19. **HSV** N (0 0 1) | refl - | vertices - | whiteCenter ~[0.0633 -0.1889 1.0000] | ref [-1 0 0] | hue0 ref  
   views: +c*|-a
20. **HSV** N (0 0 -1) | refl - | vertices - | whiteCenter ~[-0.0633 0.1889 -1.0000] | ref [1 0 0] | hue0 ref  
   views: -c*|+a
21. **HSV** N (0 0 -1) | refl - | vertices - | whiteCenter ~[-0.0633 0.1889 -1.0000] | ref [-1 0 0] | hue0 ref  
   views: -c*|-a
22. **HSV** N (0 0 1) | refl - | vertices - | whiteCenter ~[0.0633 -0.1889 1.0000] | ref [0 -1 0] | hue0 ref  
   views: +c*|-b
23. **HSV** N (0 0 -1) | refl - | vertices - | whiteCenter ~[-0.0633 0.1889 -1.0000] | ref [0 1 0] | hue0 ref  
   views: -c*|+b
24. **HSV** N (0 0 -1) | refl - | vertices - | whiteCenter ~[-0.0633 0.1889 -1.0000] | ref [0 -1 0] | hue0 ref  
   views: -c*|-b

## 211

Laue 2/m11, 2 operations, test cell [5, 6, 7] [100, 90, 90]

1. **TSL** N (0 0 1) (-1 0 0) | vertices [0 1 0] [0 -1 0] | center ~[-1.0000 0.1469 0.7253] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|-a -c*|-a +c*|+b -c*|-b  
   same colours: pymtex (max |dRGB| 6.0e-06)
2. **TSL** N (0 1 0) (1 0 0) | vertices [0 0 1] [0 0 -1] | center ~[1.0000 0.8462 0.1259] | ref [0 0 1] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a -b*|+a +b*|+c -b*|-c
3. **TSL** N (0 1 0) (-1 0 0) | vertices [0 0 1] [0 0 -1] | center ~[-1.0000 0.8462 0.1259] | ref [0 0 1] | hue0 out-of-screen (0 1 0)  
   views: +b*|-a -b*|-a +b*|-c -b*|+c
4. **TSL** N (0 0 1) (1 0 0) | vertices [0 1 0] [0 -1 0] | center ~[1.0000 0.1469 0.7253] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|+a -c*|+a +c*|-b -c*|+b

   Not an MTEX key:
  - **orix C2x TSL** N ~(0.0000 -0.1488 1.0000) (0 1 0) | vertices [1 0 0] [-1 0 0] | center ~[0.0000 1.0000 0.9930]  
    nearest MTEX key 2, max |dRGB| 1.000 on 500/500 directions > 0.01
1. **HSV** N (0 0 1) | refl ~(0.0000 1.0000 -0.2026) | vertices - | whiteCenter ~[-0.0000 1.0000 0.7399] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a +c*|+b +c*|-b -c*|+b -c*|-b  
   same colours: pymtex (max |dRGB| 4.5e-04)
2. **HSV** N ~(0.0000 1.0000 -0.2026) | refl (0 0 -1) | vertices - | whiteCenter ~[0.0000 0.9464 -1.0000] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **HSV** N ~(0.0000 1.0000 -0.2026) | refl (0 0 -1) | vertices - | whiteCenter ~[0.0000 0.9464 -1.0000] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **HSV** N (0 1 0) | refl ~(-0.0000 0.1488 -1.0000) | vertices - | whiteCenter ~[-0.0000 1.0000 -0.6953] | ref [1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a -b*|+a -b*|-a +b*|+c +b*|-c -b*|+c -b*|-c

## m11

Laue 2/m11, 2 operations, test cell [5, 6, 7] [100, 90, 90]

1. **TSL** N (0 0 1) (-1 0 0) | vertices [0 1 0] [0 -1 0] | center ~[-1.0000 0.1469 0.7253] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|-a -c*|-a +c*|+b -c*|-b  
   same colours: pymtex (max |dRGB| 6.0e-06)
2. **TSL** N (0 1 0) (1 0 0) | vertices [0 0 1] [0 0 -1] | center ~[1.0000 0.8462 0.1259] | ref [0 0 1] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a -b*|+a +b*|+c -b*|-c
3. **TSL** N (0 1 0) (-1 0 0) | vertices [0 0 1] [0 0 -1] | center ~[-1.0000 0.8462 0.1259] | ref [0 0 1] | hue0 out-of-screen (0 1 0)  
   views: +b*|-a -b*|-a +b*|-c -b*|+c
4. **TSL** N (0 0 1) (1 0 0) | vertices [0 1 0] [0 -1 0] | center ~[1.0000 0.1469 0.7253] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|+a -c*|+a +c*|-b -c*|+b

   Not an MTEX key:
  - **orix Csx TSL** N ~(0.0000 -0.1488 1.0000) (0 1 0) | vertices [1 0 0] [-1 0 0] | center ~[0.0000 1.0000 0.9930]  
    nearest MTEX key 2, max |dRGB| 1.000 on 500/500 directions > 0.01
1. **HSV** N (-1 0 0) | refl - | vertices - | whiteCenter [-1 0 0] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.4e-04)
2. **HSV** N (-1 0 0) | refl - | vertices - | whiteCenter [-1 0 0] | ref [0 1 0] | hue0 ref  
   views: +a*|+b -a*|+b
3. **HSV** N (-1 0 0) | refl - | vertices - | whiteCenter [-1 0 0] | ref [0 -1 0] | hue0 ref  
   views: +a*|-b -a*|-b
4. **HSV** N (-1 0 0) | refl - | vertices - | whiteCenter [-1 0 0] | ref [0 0 1] | hue0 ref  
   views: +a*|+c -a*|+c
5. **HSV** N (-1 0 0) | refl - | vertices - | whiteCenter [-1 0 0] | ref [0 0 -1] | hue0 ref  
   views: +a*|-c -a*|-c
6. **HSV** N (-1 0 0) | refl - | vertices - | whiteCenter [-1 0 0] | ref [1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
7. **HSV** N (-1 0 0) | refl - | vertices - | whiteCenter [-1 0 0] | ref [1 0 0] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
8. **HSV** N (-1 0 0) | refl - | vertices - | whiteCenter [-1 0 0] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 2/m11

Laue 2/m11, 4 operations, test cell [5, 6, 7] [100, 90, 90]

1. **TSL** N (0 0 1) (-1 0 0) | vertices [0 1 0] [0 -1 0] | center ~[-1.0000 0.1469 0.7253] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|-a -c*|-a +c*|+b -c*|-b  
   same colours: pymtex (max |dRGB| 6.0e-06)
2. **TSL** N (0 1 0) (1 0 0) | vertices [0 0 1] [0 0 -1] | center ~[1.0000 0.8462 0.1259] | ref [0 0 1] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a -b*|+a +b*|+c -b*|-c
3. **TSL** N (0 1 0) (-1 0 0) | vertices [0 0 1] [0 0 -1] | center ~[-1.0000 0.8462 0.1259] | ref [0 0 1] | hue0 out-of-screen (0 1 0)  
   views: +b*|-a -b*|-a +b*|-c -b*|+c
4. **TSL** N (0 0 1) (1 0 0) | vertices [0 1 0] [0 -1 0] | center ~[1.0000 0.1469 0.7253] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|+a -c*|+a +c*|-b -c*|+b

   Not an MTEX key:
  - orix: no orix point group in this setting
1. **HSV** N (0 0 1) (-1 0 0) | refl ~(-0.0000 1.0000 -0.2026) | vertices [0 1 0] [0 -1 0] | whiteCenter ~[-1.0000 0.9803 0.7253] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|-a -c*|-a +c*|+b -c*|-b  
   same colours: pymtex (max |dRGB| 4.4e-04)
2. **HSV** N (0 1 0) (1 0 0) | refl ~(0.0000 -0.1488 1.0000) | vertices [0 0 1] [0 0 -1] | whiteCenter ~[1.0000 0.8462 0.8402] | ref [1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a -b*|+a +b*|+c -b*|-c
3. **HSV** N (0 1 0) (-1 0 0) | refl ~(-0.0000 0.1488 -1.0000) | vertices [0 0 1] [0 0 -1] | whiteCenter ~[-1.0000 0.8462 -0.5883] | ref [-1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|-a -b*|-a +b*|-c -b*|+c
4. **HSV** N (0 0 1) (1 0 0) | refl ~(0.0000 -1.0000 0.2026) | vertices [0 1 0] [0 -1 0] | whiteCenter ~[1.0000 -0.6864 0.7253] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|+a -c*|+a +c*|-b -c*|+b

## 121

Laue 12/m1, 2 operations, test cell [5, 6, 7] [90, 100, 90]

1. **TSL** N (0 0 1) (0 1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2116 1.0000 0.8704] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a -c*|-a +c*|+b -c*|+b  
   same colours: pymtex (max |dRGB| 6.1e-06)
2. **TSL** N (1 0 0) (0 1 0) | vertices [0 0 1] [0 0 -1] | center ~[1.0000 0.8207 0.1240] | ref [0 0 1] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b -a*|+b +a*|-c -a*|+c
3. **TSL** N (1 0 0) (0 -1 0) | vertices [0 0 1] [0 0 -1] | center ~[1.0000 -0.8207 0.1240] | ref [0 0 1] | hue0 out-of-screen (1 0 0)  
   views: +a*|-b -a*|-b +a*|+c -a*|-c
4. **TSL** N (0 0 1) (0 -1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2116 -1.0000 0.8704] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|-a -c*|+a +c*|-b -c*|-b

   Not an MTEX key:
  - **orix C2y TSL** N ~(-0.1240 0.0000 1.0000) (0 1 0) | vertices ~[1.0000 0.0000 0.1240] ~[-1.0000 0.0000 -0.1240] | center [0 7 6]  
    nearest MTEX key 4, max |dRGB| 0.997 on 499/500 directions > 0.01
1. **HSV** N (0 0 1) | refl ~(-1.0000 -0.0000 0.2431) | vertices - | whiteCenter ~[-1.0000 -0.0000 0.8806] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a +c*|+b +c*|-b -c*|+b -c*|-b
2. **HSV** N (1 0 0) | refl ~(-0.1240 0.0000 1.0000) | vertices - | whiteCenter ~[1.0000 -0.0000 0.8275] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b -a*|+b -a*|-b +a*|+c +a*|-c -a*|+c -a*|-c

   Not an MTEX key:
  - **pymtex HSV** N (0 0 1) | vertices - | center ~[0.2431 -0.0000 1.0000]  
    nearest MTEX key 2, max |dRGB| 0.993 on 500/500 directions > 0.01

## 1m1

Laue 12/m1, 2 operations, test cell [5, 6, 7] [90, 100, 90]

1. **TSL** N (0 0 1) (0 1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2116 1.0000 0.8704] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a -c*|-a +c*|+b -c*|+b  
   same colours: pymtex (max |dRGB| 6.1e-06)
2. **TSL** N (1 0 0) (0 1 0) | vertices [0 0 1] [0 0 -1] | center ~[1.0000 0.8207 0.1240] | ref [0 0 1] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b -a*|+b +a*|-c -a*|+c
3. **TSL** N (1 0 0) (0 -1 0) | vertices [0 0 1] [0 0 -1] | center ~[1.0000 -0.8207 0.1240] | ref [0 0 1] | hue0 out-of-screen (1 0 0)  
   views: +a*|-b -a*|-b +a*|+c -a*|-c
4. **TSL** N (0 0 1) (0 -1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2116 -1.0000 0.8704] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|-a -c*|+a +c*|-b -c*|-b

   Not an MTEX key:
  - **orix Csy TSL** N ~(-0.1240 0.0000 1.0000) (0 1 0) | vertices ~[1.0000 0.0000 0.1240] ~[-1.0000 0.0000 -0.1240] | center [0 7 6]  
    nearest MTEX key 4, max |dRGB| 0.997 on 499/500 directions > 0.01
1. **HSV** N (0 1 0) | refl - | vertices - | whiteCenter [0 1 0] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.5e-04)
2. **HSV** N (0 1 0) | refl - | vertices - | whiteCenter [0 1 0] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **HSV** N (0 1 0) | refl - | vertices - | whiteCenter [0 1 0] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **HSV** N (0 1 0) | refl - | vertices - | whiteCenter [0 1 0] | ref [1 0 0] | hue0 ref  
   views: +b*|+a -b*|+a
5. **HSV** N (0 1 0) | refl - | vertices - | whiteCenter [0 1 0] | ref [-1 0 0] | hue0 ref  
   views: +b*|-a -b*|-a
6. **HSV** N (0 1 0) | refl - | vertices - | whiteCenter [0 1 0] | ref [0 0 1] | hue0 ref  
   views: +b*|+c -b*|+c
7. **HSV** N (0 1 0) | refl - | vertices - | whiteCenter [0 1 0] | ref [0 0 -1] | hue0 ref  
   views: +b*|-c -b*|-c
8. **HSV** N (0 1 0) | refl - | vertices - | whiteCenter [0 1 0] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 12/m1

Laue 12/m1, 4 operations, test cell [5, 6, 7] [90, 100, 90]

1. **TSL** N (0 0 1) (0 1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2116 1.0000 0.8704] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a -c*|-a +c*|+b -c*|+b  
   same colours: pymtex (max |dRGB| 6.1e-06)
2. **TSL** N (1 0 0) (0 1 0) | vertices [0 0 1] [0 0 -1] | center ~[1.0000 0.8207 0.1240] | ref [0 0 1] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b -a*|+b +a*|-c -a*|+c
3. **TSL** N (1 0 0) (0 -1 0) | vertices [0 0 1] [0 0 -1] | center ~[1.0000 -0.8207 0.1240] | ref [0 0 1] | hue0 out-of-screen (1 0 0)  
   views: +a*|-b -a*|-b +a*|+c -a*|-c
4. **TSL** N (0 0 1) (0 -1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2116 -1.0000 0.8704] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|-a -c*|+a +c*|-b -c*|-b

   Not an MTEX key:
  - orix: no orix point group in this setting
1. **HSV** N (0 0 1) (0 1 0) | refl ~(-1.0000 0.0000 0.2431) | vertices [-1 0 0] [1 0 0] | whiteCenter ~[-0.9884 1.0000 0.8704] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a -c*|-a +c*|+b -c*|+b
2. **HSV** N (1 0 0) (0 1 0) | refl ~(-0.1240 0.0000 1.0000) | vertices [0 0 1] [0 0 -1] | whiteCenter ~[1.0000 0.8207 0.8275] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b -a*|+b +a*|-c -a*|+c
3. **HSV** N (1 0 0) (0 -1 0) | refl ~(0.1240 -0.0000 -1.0000) | vertices [0 0 1] [0 0 -1] | whiteCenter ~[1.0000 -0.8207 -0.5794] | ref [0 -1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|-b -a*|-b +a*|+c -a*|-c
4. **HSV** N (0 0 1) (0 -1 0) | refl ~(1.0000 -0.0000 -0.2431) | vertices [-1 0 0] [1 0 0] | whiteCenter ~[1.0000 -0.7084 0.6166] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|-a -c*|+a +c*|-b -c*|-b

   Not an MTEX key:
  - **pymtex HSV** N (0 0 1) (0 1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2116 1.0000 0.8704]  
    nearest MTEX key 2, max |dRGB| 0.995 on 500/500 directions > 0.01

## 112

Laue 112/m, 2 operations, test cell [5, 6, 7] [90, 90, 100]

1. **TSL** N (0 0 1) (0 1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2084 1.0000 0.8441] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b
2. **TSL** N (0 0 -1) (0 1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2084 1.0000 -0.8441] | ref [-1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **pymtex TSL** N (-1 0 0) (0 0 1) | vertices [0 1 0] [0 -1 0] | center ~[-1.0000 -0.1447 0.7034]  
    nearest MTEX key 1, max |dRGB| 1.000 on 500/500 directions > 0.01
  - **orix C2z TSL** N (0 0 1) ~(-0.1447 1.0000 0.0000) | vertices ~[1.0000 0.1447 0.0000] ~[-1.0000 -0.1447 0.0000] | center [0 7 6]  
    nearest MTEX key 1, max |dRGB| 0.997 on 499/500 directions > 0.01
1. **HSV** N (0 1 0) | refl ~(1.0000 -0.2084 -0.0000) | vertices - | whiteCenter ~[1.0000 0.7193 -0.0000] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b
2. **HSV** N (1 0 0) | refl ~(0.1447 -1.0000 0.0000) | vertices - | whiteCenter ~[1.0000 -0.6760 -0.0000] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b -a*|+b -a*|-b +a*|+c +a*|-c -a*|+c -a*|-c
3. **HSV** N (0 1 0) | refl ~(1.0000 -0.2084 -0.0000) | vertices - | whiteCenter ~[1.0000 0.7193 -0.0000] | ref [1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a -b*|+a -b*|-a +b*|+c +b*|-c -b*|+c -b*|-c
4. **HSV** N (0 1 0) | refl ~(1.0000 -0.2084 -0.0000) | vertices - | whiteCenter ~[1.0000 0.7193 -0.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **pymtex HSV** N (-1 0 0) | vertices - | center ~[-1.0000 -0.1447 0.0000]  
    nearest MTEX key 1, max |dRGB| 0.305 on 499/500 directions > 0.01

## 11m

Laue 112/m, 2 operations, test cell [5, 6, 7] [90, 90, 100]

1. **TSL** N (0 0 1) (0 1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2084 1.0000 0.8441] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b
2. **TSL** N (0 0 -1) (0 1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2084 1.0000 -0.8441] | ref [-1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **pymtex TSL** N (-1 0 0) (0 0 1) | vertices [0 1 0] [0 -1 0] | center ~[-1.0000 -0.1447 0.7034]  
    nearest MTEX key 1, max |dRGB| 1.000 on 500/500 directions > 0.01
  - **orix Csz TSL** N (0 0 1) ~(-0.1447 1.0000 0.0000) | vertices ~[1.0000 0.1447 0.0000] ~[-1.0000 -0.1447 0.0000] | center [0 7 6]  
    nearest MTEX key 1, max |dRGB| 0.997 on 499/500 directions > 0.01
1. **HSV** N (0 0 1) | refl - | vertices - | whiteCenter [0 0 1] | ref [0 1 0] | hue0 ref  
   views: default +c*|+b -c*|+b  
   same colours: pymtex (max |dRGB| 4.5e-04)
2. **HSV** N (0 0 1) | refl - | vertices - | whiteCenter [0 0 1] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **HSV** N (0 0 1) | refl - | vertices - | whiteCenter [0 0 1] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **HSV** N (0 0 1) | refl - | vertices - | whiteCenter [0 0 1] | ref [1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
5. **HSV** N (0 0 1) | refl - | vertices - | whiteCenter [0 0 1] | ref [1 0 0] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
6. **HSV** N (0 0 1) | refl - | vertices - | whiteCenter [0 0 1] | ref [1 0 0] | hue0 ref  
   views: +c*|+a -c*|+a
7. **HSV** N (0 0 1) | refl - | vertices - | whiteCenter [0 0 1] | ref [-1 0 0] | hue0 ref  
   views: +c*|-a -c*|-a
8. **HSV** N (0 0 1) | refl - | vertices - | whiteCenter [0 0 1] | ref [0 -1 0] | hue0 ref  
   views: +c*|-b -c*|-b

## 112/m

Laue 112/m, 4 operations, test cell [5, 6, 7] [90, 90, 100]

1. **TSL** N (0 0 1) (0 1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2084 1.0000 0.8441] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b
2. **TSL** N (0 0 -1) (0 1 0) | vertices [-1 0 0] [1 0 0] | center ~[0.2084 1.0000 -0.8441] | ref [-1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **pymtex TSL** N (-1 0 0) (0 0 1) | vertices [0 1 0] [0 -1 0] | center ~[-1.0000 -0.1447 0.7034]  
    nearest MTEX key 1, max |dRGB| 1.000 on 500/500 directions > 0.01
  - **orix C2h TSL** N (0 0 1) ~(-0.1447 1.0000 0.0000) | vertices ~[1.0000 0.1447 0.0000] ~[-1.0000 -0.1447 0.0000] | center [0 7 6]  
    nearest MTEX key 1, max |dRGB| 0.997 on 499/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) | refl ~(1.0000 -0.2084 0.0000) | vertices [-1 0 0] [1 0 0] | whiteCenter ~[1.0000 0.7193 0.6072] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b
2. **HSV** N (0 0 -1) (0 1 0) | refl ~(1.0000 -0.2084 0.0000) | vertices [-1 0 0] [1 0 0] | whiteCenter ~[1.0000 0.7193 -0.6072] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **pymtex HSV** N (-1 0 0) (0 0 1) | vertices [0 1 0] [0 -1 0] | center ~[-1.0000 -0.1447 0.7034]  
    nearest MTEX key 1, max |dRGB| 0.340 on 500/500 directions > 0.01

## 222

Laue mmm, 4 operations, test cell [5, 6, 7] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [0 0 1] [0 1 0] [1 0 0] | center [42 35 30] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a  
   same colours: pymtex (max |dRGB| 3.8e-06)
2. **TSL** N (0 0 1) (0 1 0) (-1 0 0) | vertices [-1 0 0] [0 0 1] [0 1 0] | center [-42 35 30] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|+b +c*|-b -c*|+b -c*|-b

   Not an MTEX key:
  - **orix D2 TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [1 0 0] [0 0 1] [0 1 0] | center [42 35 30]  
    nearest MTEX key 1, max |dRGB| 0.008 on 0/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) | refl (1 0 0) | vertices [-1 0 0] [1 0 0] | whiteCenter [42 35 30] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (1 0 0) (0 1 0) | refl (1 0 0) | vertices [0 0 1] [0 0 -1] | whiteCenter [6 5 0] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **HSV** N (-1 0 0) (0 1 0) | refl (1 0 0) | vertices [0 0 1] [0 0 -1] | whiteCenter ~[0.0000 1.0000 0.3264] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b
4. **HSV** N (-1 0 0) (0 1 0) | refl (1 0 0) | vertices [0 0 -1] [0 0 1] | whiteCenter ~[0.0000 0.6623 1.0000] | ref [0 0 1] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+c -a*|-c
5. **HSV** N (0 0 -1) (0 1 0) | refl (1 0 0) | vertices [-1 0 0] [1 0 0] | whiteCenter [42 35 -30] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 2mm

Laue mmm, 4 operations, test cell [5, 6, 7] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [0 0 1] [0 1 0] [1 0 0] | center [42 35 30] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a  
   same colours: pymtex (max |dRGB| 3.8e-06)
2. **TSL** N (0 0 1) (0 1 0) (-1 0 0) | vertices [-1 0 0] [0 0 1] [0 1 0] | center [-42 35 30] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|+b +c*|-b -c*|+b -c*|-b

   Not an MTEX key:
  - orix: no orix point group in this setting
1. **HSV** N (0 0 1) (0 1 0) | refl - | vertices [-1 0 0] [1 0 0] | whiteCenter [0 7 6] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +b*|+a +b*|-c +c*|+a +c*|+b  
   same colours: pymtex (max |dRGB| 4.0e-04)
2. **HSV** N (0 0 1) (0 1 0) | refl - | vertices [-1 0 0] [1 0 0] | whiteCenter [0 7 6] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b
3. **HSV** N (0 0 -1) (0 1 0) | refl - | vertices [-1 0 0] [1 0 0] | whiteCenter [0 7 -6] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b
4. **HSV** N (0 0 1) (0 -1 0) | refl - | vertices [-1 0 0] [1 0 0] | whiteCenter [0 -7 6] | ref [0 0 1] | hue0 out-of-screen (1 0 0)  
   views: +a*|+c +a*|-c
5. **HSV** N (0 0 1) (0 1 0) | refl - | vertices [-1 0 0] [1 0 0] | whiteCenter [0 7 6] | ref [0 0 1] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+c -a*|-c
6. **HSV** N (0 0 1) (0 -1 0) | refl - | vertices [-1 0 0] [1 0 0] | whiteCenter [0 -7 6] | ref [-1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|-a +b*|+c +c*|-a +c*|-b
7. **HSV** N (0 0 1) (0 -1 0) | refl - | vertices [-1 0 0] [1 0 0] | whiteCenter [0 -7 6] | ref [1 0 0] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|+c -c*|+a -c*|-b
8. **HSV** N (0 0 -1) (0 -1 0) | refl - | vertices [-1 0 0] [1 0 0] | whiteCenter [0 -7 -6] | ref [-1 0 0] | hue0 out-of-screen (0 -1 0)  
   views: -b*|-a -b*|-c -c*|-a -c*|+b

## m2m

Laue mmm, 4 operations, test cell [5, 6, 7] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [0 0 1] [0 1 0] [1 0 0] | center [42 35 30] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a  
   same colours: pymtex (max |dRGB| 3.8e-06)
2. **TSL** N (0 0 1) (0 1 0) (-1 0 0) | vertices [-1 0 0] [0 0 1] [0 1 0] | center [-42 35 30] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|+b +c*|-b -c*|+b -c*|-b

   Not an MTEX key:
  - orix: no orix point group in this setting
1. **HSV** N (1 0 0) (0 0 1) | refl - | vertices [0 -1 0] [0 1 0] | whiteCenter [7 0 5] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default -a*|-b -a*|+c +c*|+a +c*|-b  
   same colours: pymtex (max |dRGB| 4.4e-04)
2. **HSV** N (1 0 0) (0 0 1) | refl - | vertices [0 -1 0] [0 1 0] | whiteCenter [7 0 5] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|+c -c*|+a -c*|+b
3. **HSV** N (1 0 0) (0 0 -1) | refl - | vertices [0 -1 0] [0 1 0] | whiteCenter [7 0 -5] | ref [0 -1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|-b +a*|-c -c*|-a -c*|-b
4. **HSV** N (1 0 0) (0 0 -1) | refl - | vertices [0 -1 0] [0 1 0] | whiteCenter [7 0 -5] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-c +c*|-a +c*|+b
5. **HSV** N (1 0 0) (0 0 -1) | refl - | vertices [0 -1 0] [0 1 0] | whiteCenter [7 0 -5] | ref [1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a
6. **HSV** N (1 0 0) (0 0 1) | refl - | vertices [0 -1 0] [0 1 0] | whiteCenter [7 0 5] | ref [1 0 0] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a
7. **HSV** N (1 0 0) (0 0 1) | refl - | vertices [0 -1 0] [0 1 0] | whiteCenter [7 0 5] | ref [0 0 1] | hue0 out-of-screen (0 1 0)  
   views: +b*|+c +b*|-c
8. **HSV** N (-1 0 0) (0 0 1) | refl - | vertices [0 -1 0] [0 1 0] | whiteCenter [-7 0 5] | ref [0 0 1] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+c -b*|-c

## mm2

Laue mmm, 4 operations, test cell [5, 6, 7] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [0 0 1] [0 1 0] [1 0 0] | center [42 35 30] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a  
   same colours: pymtex (max |dRGB| 3.8e-06)
2. **TSL** N (0 0 1) (0 1 0) (-1 0 0) | vertices [-1 0 0] [0 0 1] [0 1 0] | center [-42 35 30] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|+b +c*|-b -c*|+b -c*|-b

   Not an MTEX key:
  - **orix C2v TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [1 0 0] [0 0 1] [0 1 0] | center [42 35 30]  
    nearest MTEX key 1, max |dRGB| 0.008 on 0/500 directions > 0.01
1. **HSV** N (1 0 0) (0 1 0) | refl - | vertices [0 0 -1] [0 0 1] | whiteCenter [6 5 0] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a  
   same colours: pymtex (max |dRGB| 3.9e-04)
2. **HSV** N (1 0 0) (0 1 0) | refl - | vertices [0 0 -1] [0 0 1] | whiteCenter [6 5 0] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-c +b*|-a +b*|-c
3. **HSV** N (1 0 0) (0 -1 0) | refl - | vertices [0 0 -1] [0 0 1] | whiteCenter [6 -5 0] | ref [0 -1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|-b +a*|+c +b*|+a +b*|+c
4. **HSV** N (1 0 0) (0 1 0) | refl - | vertices [0 0 -1] [0 0 1] | whiteCenter [6 5 0] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|+c -b*|-a -b*|+c
5. **HSV** N (1 0 0) (0 -1 0) | refl - | vertices [0 0 -1] [0 0 1] | whiteCenter [6 -5 0] | ref [0 -1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|-b -a*|-c -b*|+a -b*|-c
6. **HSV** N (1 0 0) (0 -1 0) | refl - | vertices [0 0 -1] [0 0 1] | whiteCenter [6 -5 0] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a
7. **HSV** N (-1 0 0) (0 1 0) | refl - | vertices [0 0 -1] [0 0 1] | whiteCenter [-6 5 0] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|+b +c*|-b
8. **HSV** N (1 0 0) (0 1 0) | refl - | vertices [0 0 -1] [0 0 1] | whiteCenter [6 5 0] | ref [0 1 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+b -c*|-b

## mmm

Laue mmm, 8 operations, test cell [5, 6, 7] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [0 0 1] [0 1 0] [1 0 0] | center [42 35 30] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a  
   same colours: pymtex (max |dRGB| 3.8e-06)
2. **TSL** N (0 0 1) (0 1 0) (-1 0 0) | vertices [-1 0 0] [0 0 1] [0 1 0] | center [-42 35 30] | ref [-1 0 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|+b +c*|-b -c*|+b -c*|-b

   Not an MTEX key:
  - **orix D2h TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [1 0 0] [0 0 1] [0 1 0] | center [42 35 30]  
    nearest MTEX key 1, max |dRGB| 0.008 on 0/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) (1 0 0) | refl - | vertices [0 0 1] [0 1 0] [1 0 0] | whiteCenter [42 35 30] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a  
   same colours: pymtex (max |dRGB| 4.4e-04)
2. **HSV** N (0 0 1) (0 1 0) (-1 0 0) | refl - | vertices [-1 0 0] [0 0 1] [0 1 0] | whiteCenter [-42 35 30] | ref [0 1 0] | hue0 out-of-screen (0 0 1)  
   views: +c*|+b +c*|-b -c*|+b -c*|-b

## 4

Laue 4/m, 4 operations, test cell [5, 5, 9] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [0 0 1] [0 1 0] [1 0 0] | center [9 9 5] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 3.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 0 0) | vertices [0 1 0] [0 0 -1] [1 0 0] | center [9 9 -5] | ref [0 1 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix C4 TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [1 0 0] [0 0 1] [0 1 0] | center [9 9 5]  
    nearest MTEX key 1, max |dRGB| 0.008 on 0/500 directions > 0.01
1. **HSV** N (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 1] [0 0 -1] | whiteCenter ~[1.0000 0.4142 0.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.1e-04)
2. **HSV** N (1 0 0) (0 1 0) | refl (1 1 0) | vertices [0 0 1] [0 0 -1] | whiteCenter [1 1 0] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **HSV** N (-1 0 0) (0 1 0) | refl (1 1 0) | vertices [0 0 1] [0 0 -1] | whiteCenter ~[-0.4142 1.0000 0.0000] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **HSV** N (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 1] [0 0 -1] | whiteCenter ~[1.0000 0.4142 0.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## -4

Laue 4/m, 4 operations, test cell [5, 5, 9] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [0 0 1] [0 1 0] [1 0 0] | center [9 9 5] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 3.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 0 0) | vertices [0 1 0] [0 0 -1] [1 0 0] | center [9 9 -5] | ref [0 1 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix S4 TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [1 0 0] [0 0 1] [0 1 0] | center [9 9 5]  
    nearest MTEX key 1, max |dRGB| 0.008 on 0/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) | refl (1 0 0) | vertices [-1 0 0] [1 0 0] | whiteCenter [9 9 5] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.1e-04)
2. **HSV** N (1 0 0) (0 1 0) | refl (1 0 0) | vertices [0 0 1] [0 0 -1] | whiteCenter [1 1 0] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **HSV** N (-1 0 0) (0 1 0) | refl (1 0 0) | vertices [0 0 1] [0 0 -1] | whiteCenter ~[0.0000 1.0000 0.2115] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b
4. **HSV** N (-1 0 0) (0 1 0) | refl (1 0 0) | vertices [0 0 -1] [0 0 1] | whiteCenter ~[0.0000 1.0000 0.9786] | ref [0 0 1] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+c -a*|-c
5. **HSV** N (0 0 -1) (0 1 0) | refl (1 0 0) | vertices [-1 0 0] [1 0 0] | whiteCenter [9 9 -5] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 4/m

Laue 4/m, 8 operations, test cell [5, 5, 9] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [0 0 1] [0 1 0] [1 0 0] | center [9 9 5] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 3.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 0 0) | vertices [0 1 0] [0 0 -1] [1 0 0] | center [9 9 -5] | ref [0 1 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix C4h TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [1 0 0] [0 0 1] [0 1 0] | center [9 9 5]  
    nearest MTEX key 1, max |dRGB| 0.008 on 0/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 1] [0 1 0] [1 0 0] | whiteCenter ~[1.0000 0.4142 0.3254] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (0 0 -1) (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 1 0] [0 0 -1] [1 0 0] | whiteCenter ~[1.0000 0.4142 -0.3254] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 422

Laue 4/mmm, 8 operations, test cell [5, 5, 9] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [0 0 1] [1 1 0] [1 0 0] | center ~[1.0000 0.4142 0.3254] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 5.0e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -1 0) | vertices [0 0 -1] [1 1 0] [1 0 0] | center ~[1.0000 0.4142 -0.3254] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix D4 TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [1 0 0] [0 0 1] [1 1 0] | center ~[1.0000 0.4142 0.3254]  
    nearest MTEX key 1, max |dRGB| 0.006 on 0/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 1] [0 1 0] [1 0 0] | whiteCenter ~[1.0000 0.4142 0.3254] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (0 0 -1) (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 1 0] [0 0 -1] [1 0 0] | whiteCenter ~[1.0000 0.4142 -0.3254] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 4mm

Laue 4/mmm, 8 operations, test cell [5, 5, 9] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [0 0 1] [1 1 0] [1 0 0] | center ~[1.0000 0.4142 0.3254] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 5.0e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -1 0) | vertices [0 0 -1] [1 1 0] [1 0 0] | center ~[1.0000 0.4142 -0.3254] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix C4v TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [1 0 0] [0 0 1] [1 1 0] | center ~[1.0000 0.4142 0.3254]  
    nearest MTEX key 1, max |dRGB| 0.006 on 0/500 directions > 0.01
1. **HSV** N (0 1 0) (1 -1 0) | refl - | vertices [0 0 1] [0 0 -1] | whiteCenter ~[1.0000 0.4142 0.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.1e-04)
2. **HSV** N (0 1 0) (1 -1 0) | refl - | vertices [0 0 1] [0 0 -1] | whiteCenter ~[1.0000 0.4142 0.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## -42m

Laue 4/mmm, 8 operations, test cell [5, 5, 9] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [0 0 1] [1 1 0] [1 0 0] | center ~[1.0000 0.4142 0.3254] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 5.0e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -1 0) | vertices [0 0 -1] [1 1 0] [1 0 0] | center ~[1.0000 0.4142 -0.3254] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix D2d TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [1 0 0] [0 0 1] [1 1 0] | center ~[1.0000 0.4142 0.3254]  
    nearest MTEX key 1, max |dRGB| 0.006 on 0/500 directions > 0.01
1. **HSV** N (0 0 1) (1 1 0) (1 -1 0) | refl (0 1 0) | vertices [0 0 1] [1 -1 0] [1 1 0] | whiteCenter ~[1.0000 0.4142 0.3254] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (0 0 -1) (1 1 0) (1 -1 0) | refl (0 1 0) | vertices [0 0 -1] [1 -1 0] [1 1 0] | whiteCenter ~[1.0000 0.4142 -0.3254] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## -4m2

Laue 4/mmm, 8 operations, test cell [5, 5, 9] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [0 0 1] [1 1 0] [1 0 0] | center ~[1.0000 0.4142 0.3254] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 5.0e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -1 0) | vertices [0 0 -1] [1 1 0] [1 0 0] | center ~[1.0000 0.4142 -0.3254] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - orix: no orix point group in this setting
1. **HSV** N (0 0 1) (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 1] [0 1 0] [1 0 0] | whiteCenter ~[1.0000 0.4142 0.3254] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (0 0 -1) (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 1 0] [0 0 -1] [1 0 0] | whiteCenter ~[1.0000 0.4142 -0.3254] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 4/mmm

Laue 4/mmm, 16 operations, test cell [5, 5, 9] [90, 90, 90]

1. **TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [0 0 1] [1 1 0] [1 0 0] | center ~[1.0000 0.4142 0.3254] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 5.0e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -1 0) | vertices [0 0 -1] [1 1 0] [1 0 0] | center ~[1.0000 0.4142 -0.3254] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix D4h TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [1 0 0] [0 0 1] [1 1 0] | center ~[1.0000 0.4142 0.3254]  
    nearest MTEX key 1, max |dRGB| 0.006 on 0/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) (1 -1 0) | refl - | vertices [0 0 1] [1 1 0] [1 0 0] | whiteCenter ~[1.0000 0.4142 0.3254] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (0 0 -1) (0 1 0) (1 -1 0) | refl - | vertices [0 0 -1] [1 1 0] [1 0 0] | whiteCenter ~[1.0000 0.4142 -0.3254] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 3

Laue -3, 3 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [0 0 1] [0 1 0] [1 0 0] | center [9 9 5] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 0 0) | vertices [0 0 -1] [0 1 0] [1 0 0] | center [9 9 -5] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix C3 TSL** N (0 0 1) (-1 2 0) (1 1 0) | vertices [2 1 0] [0 0 1] [-1 1 0] | center ~[0.5000 1.0000 0.4811]  
    nearest MTEX key 2, max |dRGB| 0.996 on 493/500 directions > 0.01
1. **HSV** N (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 1] [0 0 -1] | whiteCenter [2 1 0] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +a*|+c -a*|+c +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.2e-04)
2. **HSV** N (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 1] [0 0 -1] | whiteCenter [2 1 0] | ref [0 1 0] | hue0 ref  
   views: +a*|+b -a*|+b +b*|+a +b*|-a +b*|+c +b*|-c
3. **HSV** N (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 1] [0 0 -1] | whiteCenter [2 1 0] | ref [0 -1 0] | hue0 ref  
   views: +a*|-b -a*|-b -b*|+a -b*|-a -b*|+c -b*|-c
4. **HSV** N (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 1] [0 0 -1] | whiteCenter [2 1 0] | ref [0 0 -1] | hue0 ref  
   views: +a*|-c -a*|-c -c*|+a -c*|-a -c*|+b -c*|-b

## -3

Laue -3, 6 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 0 0) | vertices [0 0 1] [0 1 0] [1 0 0] | center [9 9 5] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 0 0) | vertices [0 0 -1] [0 1 0] [1 0 0] | center [9 9 -5] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix S6 TSL** N (0 0 1) (-1 2 0) (1 1 0) | vertices [2 1 0] [0 0 1] [-1 1 0] | center ~[0.5000 1.0000 0.4811]  
    nearest MTEX key 2, max |dRGB| 0.996 on 493/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 1] [0 1 0] [1 0 0] | whiteCenter ~[1.0000 0.5000 0.2778] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.0e-04)
2. **HSV** N (0 0 -1) (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 -1] [0 1 0] [1 0 0] | whiteCenter ~[1.0000 0.5000 -0.2778] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 321

Laue -3m1, 6 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (1 1 0) (1 -2 0) | vertices [0 0 1] [1 -1 0] [2 1 0] | center ~[1.0000 0.0000 0.3208] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a +c*|+b +c*|-b -c*|+b -c*|-b  
   same colours: pymtex (max |dRGB| 4.7e-06)

   Not an MTEX key:
  - **orix D3x TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [1 0 0] [0 0 1] [1 1 0] | center ~[1.0000 0.5000 0.2778]  
    nearest MTEX key 1, max |dRGB| 1.000 on 496/500 directions > 0.01
1. **HSV** N (0 0 1) (1 1 0) (2 -1 0) | refl (0 1 0) (1 -1 0) | vertices [0 0 1] [1 -1 0] [1 2 0] | whiteCenter ~[1.0000 0.5000 0.2778] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (0 0 -1) (1 1 0) (2 -1 0) | refl (0 1 0) (1 -1 0) | vertices [0 0 -1] [1 -1 0] [1 2 0] | whiteCenter ~[1.0000 0.5000 -0.2778] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 312

Laue -31m, 6 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [0 0 1] [1 0 0] [1 1 0] | center ~[1.0000 0.5000 0.2778] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a +c*|+b +c*|-b -c*|+b -c*|-b  
   same colours: pymtex (max |dRGB| 3.8e-06)

   Not an MTEX key:
  - **orix D3y TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [1 0 0] [0 0 1] [1 1 0] | center ~[1.0000 0.5000 0.2778]  
    nearest MTEX key 1, max |dRGB| 0.997 on 353/500 directions > 0.01
1. **HSV** N (0 0 1) (1 1 0) (2 -1 0) | refl (1 -2 0) | vertices [0 0 1] [1 -1 0] [1 2 0] | whiteCenter ~[1.0000 -0.0000 0.3208] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.4e-04)
2. **HSV** N (0 0 -1) (1 1 0) (2 -1 0) | refl (1 -2 0) | vertices [0 0 -1] [1 -1 0] [1 2 0] | whiteCenter ~[1.0000 -0.0000 -0.3208] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 3m1

Laue -3m1, 6 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (1 1 0) (1 -2 0) | vertices [0 0 1] [1 -1 0] [2 1 0] | center ~[1.0000 0.0000 0.3208] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a +c*|+b +c*|-b -c*|+b -c*|-b  
   same colours: pymtex (max |dRGB| 4.7e-06)

   Not an MTEX key:
  - **orix C3v TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [1 0 0] [0 0 1] [1 1 0] | center ~[1.0000 0.5000 0.2778]  
    nearest MTEX key 1, max |dRGB| 1.000 on 496/500 directions > 0.01
1. **HSV** N (1 1 0) (1 -2 0) | refl - | vertices [0 0 1] [0 0 -1] | whiteCenter [1 0 0] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.2e-04)
2. **HSV** N (1 1 0) (1 -2 0) | refl - | vertices [0 0 1] [0 0 -1] | whiteCenter [1 0 0] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c +b*|+a +b*|-a +b*|+c +b*|-c
3. **HSV** N (1 1 0) (1 -2 0) | refl - | vertices [0 0 1] [0 0 -1] | whiteCenter [1 0 0] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c -b*|+a -b*|-a -b*|+c -b*|-c
4. **HSV** N (1 1 0) (1 -2 0) | refl - | vertices [0 0 1] [0 0 -1] | whiteCenter [1 0 0] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 31m

Laue -31m, 6 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [0 0 1] [1 0 0] [1 1 0] | center ~[1.0000 0.5000 0.2778] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a +c*|+b +c*|-b -c*|+b -c*|-b  
   same colours: pymtex (max |dRGB| 3.8e-06)

   Not an MTEX key:
  - **orix C3v TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [1 0 0] [0 0 1] [1 1 0] | center ~[1.0000 0.5000 0.2778]  
    nearest MTEX key 1, max |dRGB| 0.007 on 0/500 directions > 0.01
1. **HSV** N (0 1 0) (1 -1 0) | refl - | vertices [0 0 1] [0 0 -1] | whiteCenter [2 1 0] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.2e-04)
2. **HSV** N (0 1 0) (1 -1 0) | refl - | vertices [0 0 1] [0 0 -1] | whiteCenter [2 1 0] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## -3m1

Laue -3m1, 12 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (1 1 0) (1 -2 0) | vertices [0 0 1] [1 -1 0] [2 1 0] | center ~[1.0000 0.0000 0.3208] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a +c*|+b +c*|-b -c*|+b -c*|-b  
   same colours: pymtex (max |dRGB| 4.7e-06)

   Not an MTEX key:
  - **orix D3d TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [1 0 0] [0 0 1] [1 1 0] | center ~[1.0000 0.5000 0.2778]  
    nearest MTEX key 1, max |dRGB| 1.000 on 496/500 directions > 0.01
1. **HSV** N (0 0 1) (1 1 0) (1 -2 0) | refl (0 1 0) | vertices [0 0 1] [1 -1 0] [2 1 0] | whiteCenter ~[1.0000 0.2679 0.2578] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.2e-04)
2. **HSV** N (0 0 -1) (1 1 0) (1 -2 0) | refl (0 1 0) | vertices [0 0 -1] [1 -1 0] [2 1 0] | whiteCenter ~[1.0000 0.2679 -0.2578] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## -31m

Laue -31m, 12 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [0 0 1] [1 0 0] [1 1 0] | center ~[1.0000 0.5000 0.2778] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a -c*|+a -c*|-a +c*|+b +c*|-b -c*|+b -c*|-b  
   same colours: pymtex (max |dRGB| 3.8e-06)

   Not an MTEX key:
  - **orix D3d TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [1 0 0] [0 0 1] [1 1 0] | center ~[1.0000 0.5000 0.2778]  
    nearest MTEX key 1, max |dRGB| 0.007 on 0/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) (1 -1 0) | refl (1 -2 0) | vertices [0 0 1] [1 0 0] [1 1 0] | whiteCenter ~[1.0000 0.2679 0.2578] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (0 0 -1) (0 1 0) (1 -1 0) | refl (1 -2 0) | vertices [0 0 -1] [1 0 0] [1 1 0] | whiteCenter ~[1.0000 0.2679 -0.2578] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 6

Laue 6/m, 6 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [0 0 1] [1 0 0] [1 1 0] | center ~[1.0000 0.5000 0.2778] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 3.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -1 0) | vertices [0 0 -1] [1 0 0] [1 1 0] | center ~[1.0000 0.5000 -0.2778] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix C6 TSL** N (0 0 1) (-1 2 0) (2 -1 0) | vertices [2 1 0] [0 0 1] [1 2 0] | center ~[1.0000 1.0000 0.3208]  
    nearest MTEX key 2, max |dRGB| 1.000 on 493/500 directions > 0.01
1. **HSV** N (0 1 0) (1 -1 0) | refl (1 -2 0) | vertices [0 0 1] [0 0 -1] | whiteCenter ~[1.0000 0.2679 0.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 3.8e-04)
2. **HSV** N (0 1 0) (1 -1 0) | refl (1 -2 0) | vertices [0 0 1] [0 0 -1] | whiteCenter ~[1.0000 0.2679 0.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## -6

Laue 6/m, 6 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [0 0 1] [1 0 0] [1 1 0] | center ~[1.0000 0.5000 0.2778] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 3.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -1 0) | vertices [0 0 -1] [1 0 0] [1 1 0] | center ~[1.0000 0.5000 -0.2778] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix C3h TSL** N (0 0 1) (-1 2 0) (2 -1 0) | vertices [2 1 0] [0 0 1] [1 2 0] | center ~[1.0000 1.0000 0.3208]  
    nearest MTEX key 2, max |dRGB| 1.000 on 493/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 1] [0 1 0] [1 0 0] | whiteCenter ~[1.0000 0.5000 0.2778] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (0 0 -1) (0 1 0) (1 0 0) | refl (1 -1 0) | vertices [0 0 -1] [0 1 0] [1 0 0] | whiteCenter ~[1.0000 0.5000 -0.2778] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 6/m

Laue 6/m, 12 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 -1 0) | vertices [0 0 1] [1 0 0] [1 1 0] | center ~[1.0000 0.5000 0.2778] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 3.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -1 0) | vertices [0 0 -1] [1 0 0] [1 1 0] | center ~[1.0000 0.5000 -0.2778] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix C6h TSL** N (0 0 1) (-1 2 0) (2 -1 0) | vertices [2 1 0] [0 0 1] [1 2 0] | center ~[1.0000 1.0000 0.3208]  
    nearest MTEX key 2, max |dRGB| 1.000 on 493/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) (1 -1 0) | refl (1 -2 0) | vertices [0 0 1] [1 0 0] [1 1 0] | whiteCenter ~[1.0000 0.2679 0.2578] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.2e-04)
2. **HSV** N (0 0 -1) (0 1 0) (1 -1 0) | refl (1 -2 0) | vertices [0 0 -1] [1 0 0] [1 1 0] | whiteCenter ~[1.0000 0.2679 -0.2578] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 622

Laue 6/mmm, 12 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 -2 0) | vertices [0 0 1] [1 0 0] [2 1 0] | center ~[1.0000 0.2679 0.2578] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 6.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -2 0) | vertices [0 0 -1] [1 0 0] [2 1 0] | center ~[1.0000 0.2679 -0.2578] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix D6 TSL** N (0 0 1) (-1 2 0) (1 -1 0) | vertices [2 1 0] [0 0 1] [1 1 0] | center ~[1.0000 0.7321 0.2578]  
    nearest MTEX key 2, max |dRGB| 0.007 on 0/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) (1 -1 0) | refl (1 -2 0) | vertices [0 0 1] [1 0 0] [1 1 0] | whiteCenter ~[1.0000 0.2679 0.2578] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (0 0 -1) (0 1 0) (1 -1 0) | refl (1 -2 0) | vertices [0 0 -1] [1 0 0] [1 1 0] | whiteCenter ~[1.0000 0.2679 -0.2578] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 6mm

Laue 6/mmm, 12 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 -2 0) | vertices [0 0 1] [1 0 0] [2 1 0] | center ~[1.0000 0.2679 0.2578] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 6.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -2 0) | vertices [0 0 -1] [1 0 0] [2 1 0] | center ~[1.0000 0.2679 -0.2578] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix C6v TSL** N (0 0 1) (-1 2 0) (1 -1 0) | vertices [2 1 0] [0 0 1] [1 1 0] | center ~[1.0000 0.7321 0.2578]  
    nearest MTEX key 2, max |dRGB| 0.007 on 0/500 directions > 0.01
1. **HSV** N (0 1 0) (1 -2 0) | refl - | vertices [0 0 1] [0 0 -1] | whiteCenter ~[1.0000 0.2679 0.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 3.8e-04)
2. **HSV** N (0 1 0) (1 -2 0) | refl - | vertices [0 0 1] [0 0 -1] | whiteCenter ~[1.0000 0.2679 0.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## -6m2

Laue 6/mmm, 12 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 -2 0) | vertices [0 0 1] [1 0 0] [2 1 0] | center ~[1.0000 0.2679 0.2578] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 6.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -2 0) | vertices [0 0 -1] [1 0 0] [2 1 0] | center ~[1.0000 0.2679 -0.2578] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix D3h TSL** N (0 0 1) (-1 2 0) (1 -1 0) | vertices [2 1 0] [0 0 1] [1 1 0] | center ~[1.0000 0.7321 0.2578]  
    nearest MTEX key 2, max |dRGB| 0.007 on 0/500 directions > 0.01
1. **HSV** N (0 0 1) (1 1 0) (1 -2 0) | refl - | vertices [0 0 1] [1 -1 0] [2 1 0] | whiteCenter ~[1.0000 0.0000 0.3208] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.4e-04)
2. **HSV** N (0 0 -1) (1 1 0) (1 -2 0) | refl - | vertices [0 0 -1] [1 -1 0] [2 1 0] | whiteCenter ~[1.0000 0.0000 -0.3208] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## -62m

Laue 6/mmm, 12 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 -2 0) | vertices [0 0 1] [1 0 0] [2 1 0] | center ~[1.0000 0.2679 0.2578] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 6.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -2 0) | vertices [0 0 -1] [1 0 0] [2 1 0] | center ~[1.0000 0.2679 -0.2578] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix D3h TSL** N (0 0 1) (-1 2 0) (1 -1 0) | vertices [2 1 0] [0 0 1] [1 1 0] | center ~[1.0000 0.7321 0.2578]  
    nearest MTEX key 2, max |dRGB| 0.007 on 0/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) (1 -1 0) | refl - | vertices [0 0 1] [1 0 0] [1 1 0] | whiteCenter ~[1.0000 0.5000 0.2778] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (0 0 -1) (0 1 0) (1 -1 0) | refl - | vertices [0 0 -1] [1 0 0] [1 1 0] | whiteCenter ~[1.0000 0.5000 -0.2778] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 6/mmm

Laue 6/mmm, 24 operations, test cell [5, 5, 9] [90, 90, 120]

1. **TSL** N (0 0 1) (0 1 0) (1 -2 0) | vertices [0 0 1] [1 0 0] [2 1 0] | center ~[1.0000 0.2679 0.2578] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 6.8e-06)
2. **TSL** N (0 0 -1) (0 1 0) (1 -2 0) | vertices [0 0 -1] [1 0 0] [2 1 0] | center ~[1.0000 0.2679 -0.2578] | ref [0 0 -1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix D6h TSL** N (0 0 1) (-1 2 0) (1 -1 0) | vertices [2 1 0] [0 0 1] [1 1 0] | center ~[1.0000 0.7321 0.2578]  
    nearest MTEX key 2, max |dRGB| 0.007 on 0/500 directions > 0.01
1. **HSV** N (0 0 1) (0 1 0) (1 -2 0) | refl - | vertices [0 0 1] [1 0 0] [2 1 0] | whiteCenter ~[1.0000 0.2679 0.2578] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (0 0 -1) (0 1 0) (1 -2 0) | refl - | vertices [0 0 -1] [1 0 0] [2 1 0] | whiteCenter ~[1.0000 0.2679 -0.2578] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 23

Laue m-3, 12 operations, test cell [5, 5, 5] [90, 90, 90]

1. **TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | center ~[0.4014 0.4004 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b
2. **TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | center ~[0.4014 0.4004 1.0000] | ref [0 0 1] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | center ~[0.4014 0.4004 1.0000] | ref [0 0 1] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | center ~[0.4014 0.4004 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
5. **TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | center ~[0.4014 0.4004 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
6. **TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | center ~[0.4014 0.4004 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **pymtex TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [1 1 1] [0 1 1] [1 0 1] [0 0 1] | center ~[0.4024 0.3995 1.0000]  
    nearest MTEX key 1, max |dRGB| 0.003 on 0/500 directions > 0.01
  - **orix T TSL** N (1 0 0) (0 -1 1) (-1 0 1) (0 1 0) (0 0 1) | vertices [0 1 1] [1 1 1] [1 0 1] [0 0 1] | center ~[0.4023 0.4002 1.0000]  
    nearest MTEX key 1, max |dRGB| 0.002 on 0/500 directions > 0.01
1. **HSV** N (1 1 0) (1 -1 0) (0 1 1) (0 -1 1) | refl (-1 0 1) | vertices [0 0 1] [1 -1 1] [1 1 1] [1 0 0] | whiteCenter ~[0.5359 0.0000 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default -a*|+b -a*|-b -a*|+c -a*|-c +c*|+a +c*|-a +c*|+b +c*|-b
2. **HSV** N (1 1 0) (1 -1 0) (0 1 1) (0 -1 1) | refl (-1 0 1) | vertices [0 0 1] [1 -1 1] [1 1 1] [1 0 0] | whiteCenter ~[0.5359 0.0000 1.0000] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c -c*|+a -c*|-a -c*|+b -c*|-b
3. **HSV** N (1 1 0) (1 -1 0) (0 1 1) (0 -1 1) | refl (-1 0 1) | vertices [0 0 1] [1 -1 1] [1 1 1] [1 0 0] | whiteCenter ~[0.5359 0.0000 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
4. **HSV** N (1 1 0) (1 -1 0) (0 1 1) (0 -1 1) | refl (-1 0 1) | vertices [0 0 1] [1 -1 1] [1 1 1] [1 0 0] | whiteCenter ~[0.5359 0.0000 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c

   Not an MTEX key:
  - **pymtex HSV** N (1 1 0) (1 -1 0) (0 1 1) (0 -1 1) | vertices [0 0 1] [1 -1 1] [1 1 1] [1 0 0] | center ~[1.0000 -0.0000 0.9983]  
    nearest MTEX key 2, max |dRGB| 0.998 on 500/500 directions > 0.01

## m-3

Laue m-3, 24 operations, test cell [5, 5, 5] [90, 90, 90]

1. **TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | center ~[0.4014 0.4004 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b
2. **TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | center ~[0.4014 0.4004 1.0000] | ref [0 0 1] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | center ~[0.4014 0.4004 1.0000] | ref [0 0 1] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | center ~[0.4014 0.4004 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
5. **TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | center ~[0.4014 0.4004 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
6. **TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | center ~[0.4014 0.4004 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **pymtex TSL** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [1 1 1] [0 1 1] [1 0 1] [0 0 1] | center ~[0.4024 0.3995 1.0000]  
    nearest MTEX key 1, max |dRGB| 0.003 on 0/500 directions > 0.01
  - **orix Th TSL** N (1 0 0) (0 -1 1) (-1 0 1) (0 1 0) (0 0 1) | vertices [0 1 1] [1 1 1] [1 0 1] [0 0 1] | center ~[0.4023 0.4002 1.0000]  
    nearest MTEX key 1, max |dRGB| 0.002 on 0/500 directions > 0.01
1. **HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | refl (1 -1 0) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.2e-04)
2. **HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | refl (1 -1 0) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | refl (1 -1 0) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | refl (1 -1 0) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
5. **HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | refl (1 -1 0) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
6. **HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | refl (1 -1 0) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

## 432

Laue m-3m, 24 operations, test cell [5, 5, 5] [90, 90, 90]

1. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 1.0e-05)
2. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
5. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
6. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix O TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [1 1 1] [1 0 1] [0 0 1] | center ~[0.5623 0.2527 1.0000]  
    nearest MTEX key 1, max |dRGB| 0.009 on 0/500 directions > 0.01
1. **HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | refl (1 -1 0) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b
2. **HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | refl (1 -1 0) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | refl (1 -1 0) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | refl (1 -1 0) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
5. **HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | refl (1 -1 0) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
6. **HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | refl (1 -1 0) | vertices [0 0 1] [0 1 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **pymtex HSV** N (0 -1 1) (-1 0 1) (1 0 0) (0 1 0) (0 0 1) | vertices [1 1 1] [0 1 1] [1 0 1] [0 0 1] | center ~[0.4024 0.3995 1.0000]  
    nearest MTEX key 1, max |dRGB| 0.004 on 0/500 directions > 0.01

## -43m

Laue m-3m, 24 operations, test cell [5, 5, 5] [90, 90, 90]

1. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 1.0e-05)
2. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
5. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
6. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix Td TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [1 1 1] [1 0 1] [0 0 1] | center ~[0.5623 0.2527 1.0000]  
    nearest MTEX key 1, max |dRGB| 0.009 on 0/500 directions > 0.01
1. **HSV** N (1 -1 0) (1 1 0) (-1 0 1) | refl - | vertices [0 0 1] [1 -1 1] [1 1 1] | whiteCenter ~[0.5359 0.0000 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default -a*|+b -a*|-b -a*|+c -a*|-c +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.3e-04)
2. **HSV** N (1 -1 0) (1 1 0) (-1 0 1) | refl - | vertices [0 0 1] [1 -1 1] [1 1 1] | whiteCenter ~[0.5359 0.0000 1.0000] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c -c*|+a -c*|-a -c*|+b -c*|-b
3. **HSV** N (1 -1 0) (1 1 0) (-1 0 1) | refl - | vertices [0 0 1] [1 -1 1] [1 1 1] | whiteCenter ~[0.5359 0.0000 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
4. **HSV** N (1 -1 0) (1 1 0) (-1 0 1) | refl - | vertices [0 0 1] [1 -1 1] [1 1 1] | whiteCenter ~[0.5359 0.0000 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c

## m-3m

Laue m-3m, 48 operations, test cell [5, 5, 5] [90, 90, 90]

1. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 1.0e-05)
2. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
5. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
6. **TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [0 0 1] [1 1 1] [1 0 1] | center ~[0.5623 0.2527 1.0000] | ref [0 0 1] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

   Not an MTEX key:
  - **orix Oh TSL** N (1 -1 0) (-1 0 1) (0 1 0) | vertices [1 1 1] [1 0 1] [0 0 1] | center ~[0.5623 0.2527 1.0000]  
    nearest MTEX key 1, max |dRGB| 0.009 on 0/500 directions > 0.01
1. **HSV** N (1 -1 0) (-1 0 1) (0 1 0) | refl - | vertices [0 0 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 1)  
   views: default +c*|+a +c*|-a +c*|+b +c*|-b  
   same colours: pymtex (max |dRGB| 4.2e-04)
2. **HSV** N (1 -1 0) (-1 0 1) (0 1 0) | refl - | vertices [0 0 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [0 1 0] | hue0 out-of-screen (1 0 0)  
   views: +a*|+b +a*|-b +a*|+c +a*|-c
3. **HSV** N (1 -1 0) (-1 0 1) (0 1 0) | refl - | vertices [0 0 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [0 1 0] | hue0 out-of-screen (-1 0 0)  
   views: -a*|+b -a*|-b -a*|+c -a*|-c
4. **HSV** N (1 -1 0) (-1 0 1) (0 1 0) | refl - | vertices [0 0 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 1 0)  
   views: +b*|+a +b*|-a +b*|+c +b*|-c
5. **HSV** N (1 -1 0) (-1 0 1) (0 1 0) | refl - | vertices [0 0 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 -1 0)  
   views: -b*|+a -b*|-a -b*|+c -b*|-c
6. **HSV** N (1 -1 0) (-1 0 1) (0 1 0) | refl - | vertices [0 0 1] [1 1 1] [1 0 1] | whiteCenter ~[0.5623 0.2527 1.0000] | ref [1 0 0] | hue0 out-of-screen (0 0 -1)  
   views: -c*|+a -c*|-a -c*|+b -c*|-b

