#define        x              %ymm0
#define        r              %ymm1
#define        z              %ymm2
#define        y              %ymm3
#define        e              %ymm4
#define        a              %ymm5
#define        q              %ymm6
#define        c              %ymm7
#define        neg            %ymm8
#define        rng            %ymm9
#define        i              %ymm10


               .global        erfc_simd


               .text

erfc_simd:     vcmpps         $1, ZERO, x, neg
               vcvtdq2ps      neg, neg
               vmulps         NONE, neg, neg
	           vcmpps         $1, C0, x, rng
               vcvtdq2ps      rng, rng
               vmulps         NONE, rng, rng
               vpand          ABS, x, x
               vminps         C0, x, x
	           vaddps         ONE, x, q
               vmulps         q, q, q
	           vpand          I0, q, i
	           vpslld         $23, i, i
	           vsubps         I127, i, i
	           vmulps         x, x, y
	           vmovaps        y, z
	           vfmsub231ps    x, x, z
	           sub            $32, %rsp
	           vmovups        %ymm0, (%rsp)
	           vmovaps        y, %ymm0
	           vmulps         NONE, %ymm0, %ymm0
	           call           exp_simd
	           vmovaps        %ymm0, e
	           vmovups        (%rsp), %ymm0
	           add            $32, %rsp
	           vmovaps        ONE, %ymm11
	           vsubps         z, %ymm11, %ymm11
	           vmulps         e, %ymm11, e
	           lea            X0, %rax
	           vmovaps        NONE, %ymm11
	           vgatherdps     %ymm11, (%rax, i, 4), a
	           vsubps         a, x, x
	           lea            COEFF, %rax
	           vmovaps        NONE, %ymm14
	           vmovaps        %ymm14, %ymm11
	           vmovaps        %ymm14, %ymm12
	           vmovaps        %ymm14, %ymm13
	           vgatherdps     %ymm11, 252(%rax, i, 4), y
	           vgatherdps     %ymm12, 224(%rax, i, 4), %ymm15
	           vfmadd213ps    %ymm15, x, y
	           vgatherdps     %ymm13, 196(%rax, i, 4), %ymm15
	           vfmadd213ps    %ymm15, x, y
	           vmovaps        %ymm14, %ymm11
	           vmovaps        %ymm14, %ymm12
	           vmovaps        %ymm14, %ymm13
	           vgatherdps     %ymm11, 168(%rax, i, 4), %ymm15
	           vfmadd213ps    %ymm15, x, y
	           vgatherdps     %ymm12, 140(%rax, i, 4), %ymm15
	           vfmadd213ps    %ymm15, x, y
	           vgatherdps     %ymm13, 112(%rax, i, 4), %ymm15
	           vfmadd213ps    %ymm15, x, y
               vmovaps        %ymm14, %ymm11
	           vmovaps        %ymm14, %ymm12
	           vmovaps        %ymm14, %ymm13
	           vgatherdps     %ymm11, 84(%rax, i, 4), %ymm15
	           vfmadd213ps    %ymm15, x, y
	           vgatherdps     %ymm12, 56(%rax, i, 4), %ymm15
	           vfmadd213ps    %ymm15, x, y
	           vgatherdps     %ymm13, 28(%rax, i, 4), %ymm15
	           vfmadd213ps    %ymm15, x, y
	           vgatherdps     %ymm14, (%rax, i, 4), %ymm15
	           vfmadd213ps    %ymm15, x, y
	           vmulps         y, e, y
	           vmulps         y, rng, x
               ret

               .align         32
ZERO:          .float         0.0
               .float         0.0
               .float         0.0
               .float         0.0
               .float         0.0
               .float         0.0
               .float         0.0
               .float         0.0
I127:          .long          127
               .long          127
               .long          127
               .long          127
               .long          127
               .long          127
               .long          127
               .long          127
I0:            .long          0x7f800000
               .long          0x7f800000
               .long          0x7f800000
               .long          0x7f800000
               .long          0x7f800000
               .long          0x7f800000
               .long          0x7f800000
               .long          0x7f800000
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
X0:            .float         2.071067812e-01
               .float         7.071067812e-01
               .float         1.414213562e+00
               .float         2.414213562e+00
               .float         3.828427125e+00
               .float         5.828427125e+00
               .float         8.050000000e+00
C0:            .float         9.100000000
               .float         9.100000000
               .float         9.100000000
               .float         9.100000000
               .float         9.100000000
               .float         9.100000000
               .float         9.100000000
               .float         9.100000000
COEFF:         .float         8.033327565e-01
               .float         5.231565837e-01
               .float         3.362040024e-01
               .float         2.173730733e-01
               .float         1.427841750e-01
               .float         9.543357554e-02
               .float         6.955695769e-02
               .float         -7.956278443e-01
               .float         -3.885240311e-01
               .float         -1.774506471e-01
               .float         -7.880912357e-02
               .float         -3.510154995e-02
               .float         -1.592388652e-02
               .float         -8.512148338e-03
               .float         6.385528347e-01
               .float         2.484286067e-01
               .float         8.525089063e-02
               .float         2.711101839e-02
               .float         8.400449070e-03
               .float         2.622363431e-03
               .float         1.034163564e-03
               .float         -4.422528147e-01
               .float         -1.419056525e-01
               .float         -3.792512093e-02
               .float         -8.904890183e-03
               .float         -1.960695244e-03
               .float         -4.264215788e-04
               .float         -1.247544300e-04
               .float         2.734796389e-01
               .float         7.404307875e-02
               .float         1.580833512e-02
               .float         2.806355870e-03
               .float         4.470351082e-04
               .float         6.849816711e-05
               .float         1.494520147e-05
               .float         -1.542453308e-01
               .float         -3.581971576e-02
               .float         -6.227503601e-03
               .float         -8.518991127e-04
               .float         -9.970156387e-05
               .float         -1.087400143e-05
               .float         -1.778223273e-06
               .float         8.051146163e-02
               .float         1.623823828e-02
               .float         2.333771690e-03
               .float         2.498964927e-04
               .float         2.177831222e-05
               .float         1.706614075e-06
               .float         2.101680385e-07
               .float         -3.930596032e-02
               .float         -6.953584960e-03
               .float         -8.363005787e-04
               .float         -7.102731735e-05
               .float         -4.664252179e-06
               .float         -2.648930459e-07
               .float         -2.467730377e-08
               .float         1.809273268e-02
               .float         2.830327800e-03
               .float         2.877660174e-04
               .float         1.960534497e-05
               .float         9.803906664e-07
               .float         4.067606547e-08
               .float         2.878935799e-09
               .float         -7.901962821e-03
               .float         -1.100497996e-03
               .float         -9.540843869e-05
               .float         -5.265739474e-06
               .float         -2.024217687e-07
               .float         -6.181236117e-09
               .float         -3.337490199e-10



