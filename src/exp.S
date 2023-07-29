#define        x              %ymm0
#define        x0             %ymm1
#define        y              %ymm2


               .global        exp_simd


               .text

exp_simd:      vmulps         C8, x, x0
               vroundps       $0, x0, x0
	           vfnmadd231ps   C9, x0, x
	           vmovaps        C0, y
               vfmadd213ps    C1, x, y
               vfmadd213ps    C2, x, y
               vfmadd213ps    C3, x, y
               vfmadd213ps    C4, x, y
               vfmadd213ps    C5, x, y
               vfmadd213ps    C6, x, y
               vfmadd213ps    C7, x, y
               vfmadd213ps    C7, x, y
               vcvtps2dq      x0, x0
               vpaddd         I127, x0, x0
               vpslld         $23, x0, x0
               vmulps         x0, y, x
               ret

               .align         32
I127:          .long          127
               .long          127
               .long          127
               .long          127
               .long          127
               .long          127
               .long          127
               .long          127
C0:            .float         2.480158730e-05
               .float         2.480158730e-05
               .float         2.480158730e-05
               .float         2.480158730e-05
               .float         2.480158730e-05
               .float         2.480158730e-05
               .float         2.480158730e-05
               .float         2.480158730e-05
C1:            .float         1.984126984e-04
               .float         1.984126984e-04
               .float         1.984126984e-04
               .float         1.984126984e-04
               .float         1.984126984e-04
               .float         1.984126984e-04
               .float         1.984126984e-04
               .float         1.984126984e-04
C2:            .float         1.388888889e-03
               .float         1.388888889e-03
               .float         1.388888889e-03
               .float         1.388888889e-03
               .float         1.388888889e-03
               .float         1.388888889e-03
               .float         1.388888889e-03
               .float         1.388888889e-03
C3:            .float         8.333333333e-03
               .float         8.333333333e-03
               .float         8.333333333e-03
               .float         8.333333333e-03
               .float         8.333333333e-03
               .float         8.333333333e-03
               .float         8.333333333e-03
               .float         8.333333333e-03
C4:            .float         4.166666667e-02
               .float         4.166666667e-02
               .float         4.166666667e-02
               .float         4.166666667e-02
               .float         4.166666667e-02
               .float         4.166666667e-02
               .float         4.166666667e-02
               .float         4.166666667e-02
C5:            .float         1.666666667e-01
               .float         1.666666667e-01
               .float         1.666666667e-01
               .float         1.666666667e-01
               .float         1.666666667e-01
               .float         1.666666667e-01
               .float         1.666666667e-01
               .float         1.666666667e-01
C6:            .float         5.000000000e-01
               .float         5.000000000e-01
               .float         5.000000000e-01
               .float         5.000000000e-01
               .float         5.000000000e-01
               .float         5.000000000e-01
               .float         5.000000000e-01
               .float         5.000000000e-01
               .float         5.000000000e-01
C7:            .float         1.0
               .float         1.0
               .float         1.0
               .float         1.0
               .float         1.0
               .float         1.0
               .float         1.0
               .float         1.0
C8:            .float         1.442695041
               .float         1.442695041
               .float         1.442695041
               .float         1.442695041
               .float         1.442695041
               .float         1.442695041
               .float         1.442695041
               .float         1.442695041
C9:            .float         6.931471806e-01
               .float         6.931471806e-01
               .float         6.931471806e-01
               .float         6.931471806e-01
               .float         6.931471806e-01
               .float         6.931471806e-01
               .float         6.931471806e-01
               .float         6.931471806e-01
