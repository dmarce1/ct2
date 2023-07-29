#define        x              %ymm0
#define        x0             %ymm1
#define        s              %ymm2
#define        x2             %ymm3
#define        y              %ymm4


               .global        sin_simd
               .global        cos_simd


               .text

sin_simd:      vmulpd         C0, x, x0
	           vroundps       $0, x0, x0
	           vmulps         HALF, x0, s
	           vroundps       $1, s, s
	           vmulps         TWO, s, s
	           vsubps         s, x0, s
	           vmulps         NTWO, s, s
	           vaddps         ONE, s, s
	           vfnmadd231ps   PI, x0, x
	           vfnmadd231ps   C1, x0, x
	           vmulps         s, x, x
	           vmulps         x, x, x2
	           vmovaps        C2, y
	           vfmadd213ps    C3, x2, y
	           vfmadd213ps    C4, x2, y
	           vfmadd213ps    C5, x2, y
	           vfmadd213ps    C6, x2, y
	           vfmadd213ps    ONE, x2, y
	           vmulpd         x, y, x
	           ret

cos_simd:      vpand          ABS, x, x
               vmulps         C0, x, x0
               vroundps       $1, x0, x0
	           vmulps         HALF, x0, s
	           vroundps       $1, s, s
	           vmulps         TWO, s, s
	           vsubps         s, x0, s
	           vmulps         TWO, s, s
	           vsubps         ONE, s, s
	           vaddps         HALF, x0, x0
	           vfnmadd231ps   PI, x0, x
	           vfnmadd231ps   C1, x0, x
	           vmulps         s, x, x
	           vmulps         x, x, x2
	           vmovaps        C2, y
	           vfmadd213ps    C3, x2, y
	           vfmadd213ps    C4, x2, y
	           vfmadd213ps    C5, x2, y
	           vfmadd213ps    C6, x2, y
	           vfmadd213ps    ONE, x2, y
	           vmulpd         x, y, x
	           ret

               .align         32
PI:            .float         3.141592741
               .float         3.141592741
               .float         3.141592741
               .float         3.141592741
               .float         3.141592741
               .float         3.141592741
               .float         3.141592741
               .float         3.141592741
ABS:           .long          0x7fffffff
               .long          0x7fffffff
               .long          0x7fffffff
               .long          0x7fffffff
               .long          0x7fffffff
               .long          0x7fffffff
               .long          0x7fffffff
               .long          0x7fffffff
NONE:          .float         -1.0
               .float         -1.0
               .float         -1.0
               .float         -1.0
               .float         -1.0
               .float         -1.0
               .float         -1.0
               .float         -1.0
ONE:           .float         1.0
               .float         1.0
               .float         1.0
               .float         1.0
               .float         1.0
               .float         1.0
               .float         1.0
               .float         1.0
NTWO:          .float         -2.0
               .float         -2.0
               .float         -2.0
               .float         -2.0
               .float         -2.0
               .float         -2.0
               .float         -2.0
               .float         -2.0
TWO:           .float         2.0
               .float         2.0
               .float         2.0
               .float         2.0
               .float         2.0
               .float         2.0
               .float         2.0
               .float         2.0
HALF:          .float         0.5
               .float         0.5
               .float         0.5
               .float         0.5
               .float         0.5
               .float         0.5
               .float         0.5
               .float         0.5
C0:            .float         3.183098862e-01
               .float         3.183098862e-01
               .float         3.183098862e-01
               .float         3.183098862e-01
               .float         3.183098862e-01
               .float         3.183098862e-01
               .float         3.183098862e-01
               .float         3.183098862e-01
C1:            .float         -8.742277657e-08
               .float         -8.742277657e-08
               .float         -8.742277657e-08
               .float         -8.742277657e-08
               .float         -8.742277657e-08
               .float         -8.742277657e-08
               .float         -8.742277657e-08
               .float         -8.742277657e-08
C2:            .float         -2.505210794e-08
               .float         -2.505210794e-08
               .float         -2.505210794e-08
               .float         -2.505210794e-08
               .float         -2.505210794e-08
               .float         -2.505210794e-08
               .float         -2.505210794e-08
               .float         -2.505210794e-08
C3:            .float         2.755731884e-06
               .float         2.755731884e-06
               .float         2.755731884e-06
               .float         2.755731884e-06
               .float         2.755731884e-06
               .float         2.755731884e-06
               .float         2.755731884e-06
               .float         2.755731884e-06
C4:            .float         -1.984127011e-04
               .float         -1.984127011e-04
               .float         -1.984127011e-04
               .float         -1.984127011e-04
               .float         -1.984127011e-04
               .float         -1.984127011e-04
               .float         -1.984127011e-04
               .float         -1.984127011e-04
C5:            .float         8.333333768e-03
               .float         8.333333768e-03
               .float         8.333333768e-03
               .float         8.333333768e-03
               .float         8.333333768e-03
               .float         8.333333768e-03
               .float         8.333333768e-03
               .float         8.333333768e-03
C6:            .float         -1.666666716e-01
               .float         -1.666666716e-01
               .float         -1.666666716e-01
               .float         -1.666666716e-01
               .float         -1.666666716e-01
               .float         -1.666666716e-01
               .float         -1.666666716e-01
               .float         -1.666666716e-01
