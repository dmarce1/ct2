#define        x              -32(%rbp)
#define        y              -64(%rbp)
#define        z              -96(%rbp)
#define        cosphi         -128(%rbp)
#define        sinphi         -160(%rbp)
#define        R              -192(%rbp)
#define        Lsrc           -224(%rbp)
#define        Ldst           -232(%rbp)
#define        Lwrk           -240(%rbp)

               .global        L2P

               .text

L2P:           push           %rbp
               mov            %rsp, %rbp
               sub            $2296, %rsp
               and            $0xffffffffffffffc0, %rsp
               mov            %rsp, Lwrk
               sub            $224, %rsp
               mov            %rsp, %r8
               sub            $224, %rsp
               mov            %rsp, %r9
               mov            %rsi, Lsrc
               mov            %rdi, Ldst
               vmovups        %ymm0, x
               vmovups        %ymm1, y
               vmovups        %ymm2, z
               vmulps         %ymm0, %ymm0, %ymm3
               vfmadd231ps    %ymm1, %ymm1, %ymm3
               vfmadd231ps    %ymm2, %ymm2, %ymm3
               vcmpps         $1, C000, %ymm3, %ymm4
               vpand          ABS, %ymm4, %ymm4
               vaddps         %ymm3, %ymm4, %ymm5
               vrsqrtps       %ymm5, %ymm6
               vmovaps        C001, %ymm7
               vsubps         %ymm4, %ymm7, %ymm8
               vdivps         %ymm6, %ymm8, %ymm9
               vmovaps        %ymm9, R
               vfmadd231ps    x, %ymm6, %ymm4
               vmulps         y, %ymm6, %ymm5
               vmulps         C002, %ymm5, %ymm5
               vmovups        %ymm4, cosphi
               vmovups        %ymm5, sinphi
               vmovups        z, %ymm0
               mov            Lsrc, %rsi
               mov            Lsrc, %rdi
               vmovaps        %ymm0, %ymm15
               vmulps         %ymm0, %ymm0, %ymm14
               vmulps         %ymm14, %ymm15, %ymm13
               vmulps         %ymm14, %ymm14, %ymm12
               vmulps         %ymm13, %ymm14, %ymm11
               vmulps         %ymm13, %ymm13, %ymm10
               vmulps         %ymm12, %ymm13, %ymm9
               vmulps         C003, %ymm14, %ymm14
               vmulps         C004, %ymm13, %ymm13
               vmulps         C005, %ymm12, %ymm12
               vmulps         C006, %ymm11, %ymm11
               vmulps         C007, %ymm10, %ymm10
               vmulps         C008, %ymm9, %ymm9
               vmovaps        0(%rsi), %ymm0
               vmovaps        64(%rsi), %ymm1
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vmovaps        192(%rsi), %ymm1
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vmovaps        384(%rsi), %ymm1
               vfmadd231ps    %ymm13, %ymm1, %ymm0
               vmovaps        640(%rsi), %ymm1
               vfmadd231ps    %ymm12, %ymm1, %ymm0
               vmovaps        960(%rsi), %ymm1
               vfmadd231ps    %ymm11, %ymm1, %ymm0
               vmovaps        1344(%rsi), %ymm1
               vfmadd231ps    %ymm10, %ymm1, %ymm0
               vmovaps        1792(%rsi), %ymm1
               vfmadd231ps    %ymm9, %ymm1, %ymm0
               vmovaps        %ymm0, 0(%rdi)
               vmovaps        96(%rsi), %ymm0
               vmovaps        64(%rsi), %ymm2
               vmovaps        32(%rsi), %ymm4
               vmovaps        224(%rsi), %ymm1
               vmovaps        192(%rsi), %ymm3
               vmovaps        160(%rsi), %ymm5
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vmovaps        416(%rsi), %ymm1
               vmovaps        384(%rsi), %ymm3
               vmovaps        352(%rsi), %ymm5
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vfmadd231ps    %ymm14, %ymm3, %ymm2
               vfmadd231ps    %ymm14, %ymm5, %ymm4
               vmovaps        672(%rsi), %ymm1
               vmovaps        640(%rsi), %ymm3
               vmovaps        608(%rsi), %ymm5
               vfmadd231ps    %ymm13, %ymm1, %ymm0
               vfmadd231ps    %ymm13, %ymm3, %ymm2
               vfmadd231ps    %ymm13, %ymm5, %ymm4
               vmovaps        992(%rsi), %ymm1
               vmovaps        960(%rsi), %ymm3
               vmovaps        928(%rsi), %ymm5
               vfmadd231ps    %ymm12, %ymm1, %ymm0
               vfmadd231ps    %ymm12, %ymm3, %ymm2
               vfmadd231ps    %ymm12, %ymm5, %ymm4
               vmovaps        1376(%rsi), %ymm1
               vmovaps        1344(%rsi), %ymm3
               vmovaps        1312(%rsi), %ymm5
               vfmadd231ps    %ymm11, %ymm1, %ymm0
               vfmadd231ps    %ymm11, %ymm3, %ymm2
               vfmadd231ps    %ymm11, %ymm5, %ymm4
               vmovaps        1824(%rsi), %ymm1
               vmovaps        1792(%rsi), %ymm3
               vmovaps        1760(%rsi), %ymm5
               vfmadd231ps    %ymm10, %ymm1, %ymm0
               vfmadd231ps    %ymm10, %ymm3, %ymm2
               vfmadd231ps    %ymm10, %ymm5, %ymm4
               vmovaps        %ymm0, 96(%rdi)
               vmovaps        %ymm2, 64(%rdi)
               vmovaps        %ymm4, 32(%rdi)
               vmovaps        256(%rsi), %ymm0
               vmovaps        224(%rsi), %ymm2
               vmovaps        192(%rsi), %ymm4
               vmovaps        160(%rsi), %ymm6
               vmovaps        448(%rsi), %ymm1
               vmovaps        416(%rsi), %ymm3
               vmovaps        384(%rsi), %ymm5
               vmovaps        352(%rsi), %ymm7
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vfmadd231ps    %ymm15, %ymm7, %ymm6
               vmovaps        704(%rsi), %ymm1
               vmovaps        672(%rsi), %ymm3
               vmovaps        640(%rsi), %ymm5
               vmovaps        608(%rsi), %ymm7
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vfmadd231ps    %ymm14, %ymm3, %ymm2
               vfmadd231ps    %ymm14, %ymm5, %ymm4
               vfmadd231ps    %ymm14, %ymm7, %ymm6
               vmovaps        1024(%rsi), %ymm1
               vmovaps        992(%rsi), %ymm3
               vmovaps        960(%rsi), %ymm5
               vmovaps        928(%rsi), %ymm7
               vfmadd231ps    %ymm13, %ymm1, %ymm0
               vfmadd231ps    %ymm13, %ymm3, %ymm2
               vfmadd231ps    %ymm13, %ymm5, %ymm4
               vfmadd231ps    %ymm13, %ymm7, %ymm6
               vmovaps        1408(%rsi), %ymm1
               vmovaps        1376(%rsi), %ymm3
               vmovaps        1344(%rsi), %ymm5
               vmovaps        1312(%rsi), %ymm7
               vfmadd231ps    %ymm12, %ymm1, %ymm0
               vfmadd231ps    %ymm12, %ymm3, %ymm2
               vfmadd231ps    %ymm12, %ymm5, %ymm4
               vfmadd231ps    %ymm12, %ymm7, %ymm6
               vmovaps        1856(%rsi), %ymm1
               vmovaps        1824(%rsi), %ymm3
               vmovaps        1792(%rsi), %ymm5
               vmovaps        1760(%rsi), %ymm7
               vfmadd231ps    %ymm11, %ymm1, %ymm0
               vfmadd231ps    %ymm11, %ymm3, %ymm2
               vfmadd231ps    %ymm11, %ymm5, %ymm4
               vfmadd231ps    %ymm11, %ymm7, %ymm6
               vmovaps        %ymm0, 256(%rdi)
               vmovaps        %ymm2, 224(%rdi)
               vmovaps        %ymm4, 192(%rdi)
               vmovaps        %ymm6, 160(%rdi)
               vmovaps        128(%rsi), %ymm0
               vmovaps        320(%rsi), %ymm1
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vmovaps        576(%rsi), %ymm1
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vmovaps        896(%rsi), %ymm1
               vfmadd231ps    %ymm13, %ymm1, %ymm0
               vmovaps        1280(%rsi), %ymm1
               vfmadd231ps    %ymm12, %ymm1, %ymm0
               vmovaps        1728(%rsi), %ymm1
               vfmadd231ps    %ymm11, %ymm1, %ymm0
               vmovaps        %ymm0, 128(%rdi)
               vmovaps        480(%rsi), %ymm0
               vmovaps        448(%rsi), %ymm2
               vmovaps        416(%rsi), %ymm4
               vmovaps        384(%rsi), %ymm6
               vmovaps        736(%rsi), %ymm1
               vmovaps        704(%rsi), %ymm3
               vmovaps        672(%rsi), %ymm5
               vmovaps        640(%rsi), %ymm7
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vfmadd231ps    %ymm15, %ymm7, %ymm6
               vmovaps        1056(%rsi), %ymm1
               vmovaps        1024(%rsi), %ymm3
               vmovaps        992(%rsi), %ymm5
               vmovaps        960(%rsi), %ymm7
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vfmadd231ps    %ymm14, %ymm3, %ymm2
               vfmadd231ps    %ymm14, %ymm5, %ymm4
               vfmadd231ps    %ymm14, %ymm7, %ymm6
               vmovaps        1440(%rsi), %ymm1
               vmovaps        1408(%rsi), %ymm3
               vmovaps        1376(%rsi), %ymm5
               vmovaps        1344(%rsi), %ymm7
               vfmadd231ps    %ymm13, %ymm1, %ymm0
               vfmadd231ps    %ymm13, %ymm3, %ymm2
               vfmadd231ps    %ymm13, %ymm5, %ymm4
               vfmadd231ps    %ymm13, %ymm7, %ymm6
               vmovaps        1888(%rsi), %ymm1
               vmovaps        1856(%rsi), %ymm3
               vmovaps        1824(%rsi), %ymm5
               vmovaps        1792(%rsi), %ymm7
               vfmadd231ps    %ymm12, %ymm1, %ymm0
               vfmadd231ps    %ymm12, %ymm3, %ymm2
               vfmadd231ps    %ymm12, %ymm5, %ymm4
               vfmadd231ps    %ymm12, %ymm7, %ymm6
               vmovaps        %ymm0, 480(%rdi)
               vmovaps        %ymm2, 448(%rdi)
               vmovaps        %ymm4, 416(%rdi)
               vmovaps        %ymm6, 384(%rdi)
               vmovaps        352(%rsi), %ymm0
               vmovaps        320(%rsi), %ymm2
               vmovaps        288(%rsi), %ymm4
               vmovaps        608(%rsi), %ymm1
               vmovaps        576(%rsi), %ymm3
               vmovaps        544(%rsi), %ymm5
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vmovaps        928(%rsi), %ymm1
               vmovaps        896(%rsi), %ymm3
               vmovaps        864(%rsi), %ymm5
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vfmadd231ps    %ymm14, %ymm3, %ymm2
               vfmadd231ps    %ymm14, %ymm5, %ymm4
               vmovaps        1312(%rsi), %ymm1
               vmovaps        1280(%rsi), %ymm3
               vmovaps        1248(%rsi), %ymm5
               vfmadd231ps    %ymm13, %ymm1, %ymm0
               vfmadd231ps    %ymm13, %ymm3, %ymm2
               vfmadd231ps    %ymm13, %ymm5, %ymm4
               vmovaps        1760(%rsi), %ymm1
               vmovaps        1728(%rsi), %ymm3
               vmovaps        1696(%rsi), %ymm5
               vfmadd231ps    %ymm12, %ymm1, %ymm0
               vfmadd231ps    %ymm12, %ymm3, %ymm2
               vfmadd231ps    %ymm12, %ymm5, %ymm4
               vmovaps        %ymm0, 352(%rdi)
               vmovaps        %ymm2, 320(%rdi)
               vmovaps        %ymm4, 288(%rdi)
               vmovaps        768(%rsi), %ymm0
               vmovaps        736(%rsi), %ymm2
               vmovaps        704(%rsi), %ymm4
               vmovaps        672(%rsi), %ymm6
               vmovaps        1088(%rsi), %ymm1
               vmovaps        1056(%rsi), %ymm3
               vmovaps        1024(%rsi), %ymm5
               vmovaps        992(%rsi), %ymm7
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vfmadd231ps    %ymm15, %ymm7, %ymm6
               vmovaps        1472(%rsi), %ymm1
               vmovaps        1440(%rsi), %ymm3
               vmovaps        1408(%rsi), %ymm5
               vmovaps        1376(%rsi), %ymm7
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vfmadd231ps    %ymm14, %ymm3, %ymm2
               vfmadd231ps    %ymm14, %ymm5, %ymm4
               vfmadd231ps    %ymm14, %ymm7, %ymm6
               vmovaps        1920(%rsi), %ymm1
               vmovaps        1888(%rsi), %ymm3
               vmovaps        1856(%rsi), %ymm5
               vmovaps        1824(%rsi), %ymm7
               vfmadd231ps    %ymm13, %ymm1, %ymm0
               vfmadd231ps    %ymm13, %ymm3, %ymm2
               vfmadd231ps    %ymm13, %ymm5, %ymm4
               vfmadd231ps    %ymm13, %ymm7, %ymm6
               vmovaps        %ymm0, 768(%rdi)
               vmovaps        %ymm2, 736(%rdi)
               vmovaps        %ymm4, 704(%rdi)
               vmovaps        %ymm6, 672(%rdi)
               vmovaps        640(%rsi), %ymm0
               vmovaps        608(%rsi), %ymm2
               vmovaps        576(%rsi), %ymm4
               vmovaps        544(%rsi), %ymm6
               vmovaps        960(%rsi), %ymm1
               vmovaps        928(%rsi), %ymm3
               vmovaps        896(%rsi), %ymm5
               vmovaps        864(%rsi), %ymm7
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vfmadd231ps    %ymm15, %ymm7, %ymm6
               vmovaps        1344(%rsi), %ymm1
               vmovaps        1312(%rsi), %ymm3
               vmovaps        1280(%rsi), %ymm5
               vmovaps        1248(%rsi), %ymm7
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vfmadd231ps    %ymm14, %ymm3, %ymm2
               vfmadd231ps    %ymm14, %ymm5, %ymm4
               vfmadd231ps    %ymm14, %ymm7, %ymm6
               vmovaps        1792(%rsi), %ymm1
               vmovaps        1760(%rsi), %ymm3
               vmovaps        1728(%rsi), %ymm5
               vmovaps        1696(%rsi), %ymm7
               vfmadd231ps    %ymm13, %ymm1, %ymm0
               vfmadd231ps    %ymm13, %ymm3, %ymm2
               vfmadd231ps    %ymm13, %ymm5, %ymm4
               vfmadd231ps    %ymm13, %ymm7, %ymm6
               vmovaps        %ymm0, 640(%rdi)
               vmovaps        %ymm2, 608(%rdi)
               vmovaps        %ymm4, 576(%rdi)
               vmovaps        %ymm6, 544(%rdi)
               vmovaps        512(%rsi), %ymm0
               vmovaps        832(%rsi), %ymm1
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vmovaps        1216(%rsi), %ymm1
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vmovaps        1664(%rsi), %ymm1
               vfmadd231ps    %ymm13, %ymm1, %ymm0
               vmovaps        %ymm0, 512(%rdi)
               vmovaps        1120(%rsi), %ymm0
               vmovaps        1088(%rsi), %ymm2
               vmovaps        1056(%rsi), %ymm4
               vmovaps        1024(%rsi), %ymm6
               vmovaps        1504(%rsi), %ymm1
               vmovaps        1472(%rsi), %ymm3
               vmovaps        1440(%rsi), %ymm5
               vmovaps        1408(%rsi), %ymm7
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vfmadd231ps    %ymm15, %ymm7, %ymm6
               vmovaps        1952(%rsi), %ymm1
               vmovaps        1920(%rsi), %ymm3
               vmovaps        1888(%rsi), %ymm5
               vmovaps        1856(%rsi), %ymm7
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vfmadd231ps    %ymm14, %ymm3, %ymm2
               vfmadd231ps    %ymm14, %ymm5, %ymm4
               vfmadd231ps    %ymm14, %ymm7, %ymm6
               vmovaps        %ymm0, 1120(%rdi)
               vmovaps        %ymm2, 1088(%rdi)
               vmovaps        %ymm4, 1056(%rdi)
               vmovaps        %ymm6, 1024(%rdi)
               vmovaps        992(%rsi), %ymm0
               vmovaps        960(%rsi), %ymm2
               vmovaps        928(%rsi), %ymm4
               vmovaps        896(%rsi), %ymm6
               vmovaps        1376(%rsi), %ymm1
               vmovaps        1344(%rsi), %ymm3
               vmovaps        1312(%rsi), %ymm5
               vmovaps        1280(%rsi), %ymm7
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vfmadd231ps    %ymm15, %ymm7, %ymm6
               vmovaps        1824(%rsi), %ymm1
               vmovaps        1792(%rsi), %ymm3
               vmovaps        1760(%rsi), %ymm5
               vmovaps        1728(%rsi), %ymm7
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vfmadd231ps    %ymm14, %ymm3, %ymm2
               vfmadd231ps    %ymm14, %ymm5, %ymm4
               vfmadd231ps    %ymm14, %ymm7, %ymm6
               vmovaps        %ymm0, 992(%rdi)
               vmovaps        %ymm2, 960(%rdi)
               vmovaps        %ymm4, 928(%rdi)
               vmovaps        %ymm6, 896(%rdi)
               vmovaps        864(%rsi), %ymm0
               vmovaps        832(%rsi), %ymm2
               vmovaps        800(%rsi), %ymm4
               vmovaps        1248(%rsi), %ymm1
               vmovaps        1216(%rsi), %ymm3
               vmovaps        1184(%rsi), %ymm5
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vmovaps        1696(%rsi), %ymm1
               vmovaps        1664(%rsi), %ymm3
               vmovaps        1632(%rsi), %ymm5
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vfmadd231ps    %ymm14, %ymm3, %ymm2
               vfmadd231ps    %ymm14, %ymm5, %ymm4
               vmovaps        %ymm0, 864(%rdi)
               vmovaps        %ymm2, 832(%rdi)
               vmovaps        %ymm4, 800(%rdi)
               vmovaps        1536(%rsi), %ymm0
               vmovaps        1504(%rsi), %ymm2
               vmovaps        1472(%rsi), %ymm4
               vmovaps        1440(%rsi), %ymm6
               vmovaps        1984(%rsi), %ymm1
               vmovaps        1952(%rsi), %ymm3
               vmovaps        1920(%rsi), %ymm5
               vmovaps        1888(%rsi), %ymm7
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vfmadd231ps    %ymm15, %ymm7, %ymm6
               vmovaps        %ymm0, 1536(%rdi)
               vmovaps        %ymm2, 1504(%rdi)
               vmovaps        %ymm4, 1472(%rdi)
               vmovaps        %ymm6, 1440(%rdi)
               vmovaps        1408(%rsi), %ymm0
               vmovaps        1376(%rsi), %ymm2
               vmovaps        1344(%rsi), %ymm4
               vmovaps        1312(%rsi), %ymm6
               vmovaps        1856(%rsi), %ymm1
               vmovaps        1824(%rsi), %ymm3
               vmovaps        1792(%rsi), %ymm5
               vmovaps        1760(%rsi), %ymm7
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vfmadd231ps    %ymm15, %ymm7, %ymm6
               vmovaps        %ymm0, 1408(%rdi)
               vmovaps        %ymm2, 1376(%rdi)
               vmovaps        %ymm4, 1344(%rdi)
               vmovaps        %ymm6, 1312(%rdi)
               vmovaps        1280(%rsi), %ymm0
               vmovaps        1248(%rsi), %ymm2
               vmovaps        1216(%rsi), %ymm4
               vmovaps        1184(%rsi), %ymm6
               vmovaps        1728(%rsi), %ymm1
               vmovaps        1696(%rsi), %ymm3
               vmovaps        1664(%rsi), %ymm5
               vmovaps        1632(%rsi), %ymm7
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vfmadd231ps    %ymm15, %ymm7, %ymm6
               vmovaps        %ymm0, 1280(%rdi)
               vmovaps        %ymm2, 1248(%rdi)
               vmovaps        %ymm4, 1216(%rdi)
               vmovaps        %ymm6, 1184(%rdi)
               vmovaps        1152(%rsi), %ymm0
               vmovaps        1600(%rsi), %ymm1
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vmovaps        %ymm0, 1152(%rdi)
               mov            Lwrk, %rsi
               mov            Ldst, %rdi
               vmovaps        %ymm0, (%r8)
               vmovaps        %ymm1, (%r9)
               vmovaps        %ymm0, %ymm4
               vmovaps        %ymm1, %ymm5
               vmulps         %ymm5, %ymm1, %ymm2
               vmulps         %ymm5, %ymm0, %ymm3
               vfnmadd231ps   %ymm4, %ymm0, %ymm2
               vfmadd231ps    %ymm4, %ymm1, %ymm3
               vmovaps        %ymm2, 32(%r8)
               vmovaps        %ymm3, 32(%r9)
               vmovaps        %ymm2, %ymm4
               vmovaps        %ymm3, %ymm5
               vmulps         %ymm5, %ymm1, %ymm2
               vmulps         %ymm5, %ymm0, %ymm3
               vfnmadd231ps   %ymm4, %ymm0, %ymm2
               vfmadd231ps    %ymm4, %ymm1, %ymm3
               vmovaps        %ymm2, 64(%r8)
               vmovaps        %ymm3, 64(%r9)
               vmovaps        %ymm2, %ymm4
               vmovaps        %ymm3, %ymm5
               vmulps         %ymm5, %ymm1, %ymm2
               vmulps         %ymm5, %ymm0, %ymm3
               vfnmadd231ps   %ymm4, %ymm0, %ymm2
               vfmadd231ps    %ymm4, %ymm1, %ymm3
               vmovaps        %ymm2, 96(%r8)
               vmovaps        %ymm3, 96(%r9)
               vmovaps        %ymm2, %ymm4
               vmovaps        %ymm3, %ymm5
               vmulps         %ymm5, %ymm1, %ymm2
               vmulps         %ymm5, %ymm0, %ymm3
               vfnmadd231ps   %ymm4, %ymm0, %ymm2
               vfmadd231ps    %ymm4, %ymm1, %ymm3
               vmovaps        %ymm2, 128(%r8)
               vmovaps        %ymm3, 128(%r9)
               vmovaps        %ymm2, %ymm4
               vmovaps        %ymm3, %ymm5
               vmulps         %ymm5, %ymm1, %ymm2
               vmulps         %ymm5, %ymm0, %ymm3
               vfnmadd231ps   %ymm4, %ymm0, %ymm2
               vfmadd231ps    %ymm4, %ymm1, %ymm3
               vmovaps        %ymm2, 160(%r8)
               vmovaps        %ymm3, 160(%r9)
               vmovaps        %ymm2, %ymm4
               vmovaps        %ymm3, %ymm5
               vmulps         %ymm5, %ymm1, %ymm2
               vmulps         %ymm5, %ymm0, %ymm3
               vfnmadd231ps   %ymm4, %ymm0, %ymm2
               vfmadd231ps    %ymm4, %ymm1, %ymm3
               vmovaps        %ymm2, 192(%r8)
               vmovaps        %ymm3, 192(%r9)
               vmovaps        0(%r8), %ymm14
               vmovaps        0(%r9), %ymm15
               vmovaps        96(%rsi), %ymm2
               vmovaps        32(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 96(%rdi)
               vmovaps        %ymm1, 32(%rdi)
               vmulps         160(%rsi), %ymm15, %ymm1
               vfmsub231ps    224(%rsi), %ymm14, %ymm1
               vmovaps        %ymm1, 224(%rdi)
               vmovaps        416(%rsi), %ymm2
               vmovaps        352(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 416(%rdi)
               vmovaps        %ymm1, 352(%rdi)
               vmulps         608(%rsi), %ymm15, %ymm1
               vfmsub231ps    672(%rsi), %ymm14, %ymm1
               vmovaps        %ymm1, 672(%rdi)
               vmovaps        992(%rsi), %ymm2
               vmovaps        928(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 992(%rdi)
               vmovaps        %ymm1, 928(%rdi)
               vmulps         1312(%rsi), %ymm15, %ymm1
               vfmsub231ps    1376(%rsi), %ymm14, %ymm1
               vmovaps        %ymm1, 1376(%rdi)
               vmovaps        1824(%rsi), %ymm2
               vmovaps        1760(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 1824(%rdi)
               vmovaps        %ymm1, 1760(%rdi)
               vmovaps        32(%r8), %ymm14
               vmovaps        32(%r9), %ymm15
               vmovaps        256(%rsi), %ymm2
               vmovaps        128(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 256(%rdi)
               vmovaps        %ymm1, 128(%rdi)
               vmulps         320(%rsi), %ymm15, %ymm1
               vfmsub231ps    448(%rsi), %ymm14, %ymm1
               vmovaps        %ymm1, 448(%rdi)
               vmovaps        704(%rsi), %ymm2
               vmovaps        576(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 704(%rdi)
               vmovaps        %ymm1, 576(%rdi)
               vmulps         896(%rsi), %ymm15, %ymm1
               vfmsub231ps    1024(%rsi), %ymm14, %ymm1
               vmovaps        %ymm1, 1024(%rdi)
               vmovaps        1408(%rsi), %ymm2
               vmovaps        1280(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 1408(%rdi)
               vmovaps        %ymm1, 1280(%rdi)
               vmulps         1728(%rsi), %ymm15, %ymm1
               vfmsub231ps    1856(%rsi), %ymm14, %ymm1
               vmovaps        %ymm1, 1856(%rdi)
               vmovaps        64(%r8), %ymm14
               vmovaps        64(%r9), %ymm15
               vmovaps        480(%rsi), %ymm2
               vmovaps        288(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 480(%rdi)
               vmovaps        %ymm1, 288(%rdi)
               vmulps         544(%rsi), %ymm15, %ymm1
               vfmsub231ps    736(%rsi), %ymm14, %ymm1
               vmovaps        %ymm1, 736(%rdi)
               vmovaps        1056(%rsi), %ymm2
               vmovaps        864(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 1056(%rdi)
               vmovaps        %ymm1, 864(%rdi)
               vmulps         1248(%rsi), %ymm15, %ymm1
               vfmsub231ps    1440(%rsi), %ymm14, %ymm1
               vmovaps        %ymm1, 1440(%rdi)
               vmovaps        1888(%rsi), %ymm2
               vmovaps        1696(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 1888(%rdi)
               vmovaps        %ymm1, 1696(%rdi)
               vmovaps        96(%r8), %ymm14
               vmovaps        96(%r9), %ymm15
               vmovaps        768(%rsi), %ymm2
               vmovaps        512(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 768(%rdi)
               vmovaps        %ymm1, 512(%rdi)
               vmulps         832(%rsi), %ymm15, %ymm1
               vfmsub231ps    1088(%rsi), %ymm14, %ymm1
               vmovaps        %ymm1, 1088(%rdi)
               vmovaps        1472(%rsi), %ymm2
               vmovaps        1216(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 1472(%rdi)
               vmovaps        %ymm1, 1216(%rdi)
               vmulps         1664(%rsi), %ymm15, %ymm1
               vfmsub231ps    1920(%rsi), %ymm14, %ymm1
               vmovaps        %ymm1, 1920(%rdi)
               vmovaps        128(%r8), %ymm14
               vmovaps        128(%r9), %ymm15
               vmovaps        1120(%rsi), %ymm2
               vmovaps        800(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 1120(%rdi)
               vmovaps        %ymm1, 800(%rdi)
               vmulps         1184(%rsi), %ymm15, %ymm1
               vfmsub231ps    1504(%rsi), %ymm14, %ymm1
               vmovaps        %ymm1, 1504(%rdi)
               vmovaps        1952(%rsi), %ymm2
               vmovaps        1632(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 1952(%rdi)
               vmovaps        %ymm1, 1632(%rdi)
               vmovaps        160(%r8), %ymm14
               vmovaps        160(%r9), %ymm15
               vmovaps        1536(%rsi), %ymm2
               vmovaps        1152(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 1536(%rdi)
               vmovaps        %ymm1, 1152(%rdi)
               vmulps         1600(%rsi), %ymm15, %ymm1
               vfmsub231ps    1984(%rsi), %ymm14, %ymm1
               vmovaps        %ymm1, 1984(%rdi)
               vmovaps        192(%r8), %ymm14
               vmovaps        192(%r9), %ymm15
               vmovaps        2016(%rsi), %ymm2
               vmovaps        1568(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmsub231ps    %ymm2, %ymm14, %ymm0
               vfmadd231ps    %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 2016(%rdi)
               vmovaps        %ymm1, 1568(%rdi)
               vmovups        cosphi, %ymm0
               vmovups        sinphi, %ymm1
               mov            Lwrk, %rdi
               mov            Ldst, %rsi
               vmovaps        1760(%rsi), %ymm1
               vfmadd231ps    C009, %ymm1, %ymm0
               vmovaps        32(%rsi), %ymm2
               vmovaps        96(%rsi), %ymm4
               vaddps         64(%rsi), %ymm6, %ymm6
               vmovaps        1792(%rsi), %ymm9
               vmulps         C010, %ymm9, %ymm8
               vmovaps        1824(%rsi), %ymm11
               vfmadd231ps    C011, %ymm11, %ymm10
               vmovaps        192(%rsi), %ymm13
               vfmadd231ps    C012, %ymm13, %ymm12
               vaddps         224(%rsi), %ymm14, %ymm14
               vmovaps        1696(%rsi), %ymm1
               vfmadd231ps    C013, %ymm1, %ymm0
               vmovaps        %ymm2, 32(%rdi)
               vmovaps        %ymm4, 64(%rdi)
               vmovaps        %ymm6, 96(%rdi)
               vmovaps        1856(%rsi), %ymm9
               vmulps         C014, %ymm9, %ymm8
               vmovaps        1888(%rsi), %ymm11
               vfmadd231ps    C015, %ymm11, %ymm10
               vmovaps        256(%rsi), %ymm13
               vfmadd231ps    C003, %ymm13, %ymm12
               vmovaps        %ymm14, 224(%rdi)
               vmovaps        1632(%rsi), %ymm1
               vfmadd231ps    C009, %ymm1, %ymm0
               vmovaps        1280(%rsi), %ymm3
               vmulps         C016, %ymm3, %ymm2
               vmovaps        1376(%rsi), %ymm5
               vfmadd231ps    C017, %ymm5, %ymm4
               vmovaps        1344(%rsi), %ymm7
               vfmadd231ps    C010, %ymm7, %ymm6
               vmovaps        1920(%rsi), %ymm9
               vmulps         C018, %ymm9, %ymm8
               vmovaps        1952(%rsi), %ymm11
               vfmadd231ps    C019, %ymm11, %ymm10
               vmovaps        128(%rsi), %ymm13
               vfmadd231ps    C003, %ymm13, %ymm12
               vmovaps        928(%rsi), %ymm15
               vmulps         C020, %ymm15, %ymm14
               vmovaps        1568(%rsi), %ymm1
               vfmadd231ps    C021, %ymm1, %ymm0
               vmovaps        1216(%rsi), %ymm3
               vmulps         C022, %ymm3, %ymm2
               vmovaps        1440(%rsi), %ymm5
               vfmadd231ps    C010, %ymm5, %ymm4
               vmovaps        1408(%rsi), %ymm7
               vfmadd231ps    C014, %ymm7, %ymm6
               vmovaps        1984(%rsi), %ymm9
               vmulps         C023, %ymm9, %ymm8
               vmovaps        2016(%rsi), %ymm11
               vfmadd231ps    C021, %ymm11, %ymm10
               vmovaps        %ymm12, 192(%rdi)
               vmovaps        864(%rsi), %ymm15
               vmulps         C018, %ymm15, %ymm14
               vmovaps        %ymm0, 1952(%rdi)
               vmovaps        1152(%rsi), %ymm3
               vmulps         C023, %ymm3, %ymm2
               vmovaps        1504(%rsi), %ymm5
               vfmadd231ps    C024, %ymm5, %ymm4
               vmovaps        1472(%rsi), %ymm7
               vfmadd231ps    C018, %ymm7, %ymm6
               vmovaps        %ymm8, 1824(%rdi)
               vmovaps        %ymm10, 1792(%rdi)
               vmovaps        960(%rsi), %ymm13
               vfmadd231ps    C025, %ymm13, %ymm12
               vmovaps        800(%rsi), %ymm15
               vmulps         C024, %ymm15, %ymm14
               vmovaps        992(%rsi), %ymm1
               vmulps         C017, %ymm1, %ymm0
               vmovaps        %ymm2, 1472(%rdi)
               vmovaps        %ymm4, 1376(%rdi)
               vmovaps        1536(%rsi), %ymm7
               vfmadd231ps    C023, %ymm7, %ymm6
               vmovaps        576(%rsi), %ymm9
               vfmadd231ps    C026, %ymm9, %ymm8
               vmovaps        672(%rsi), %ymm11
               vmulps         C027, %ymm11, %ymm10
               vmovaps        1024(%rsi), %ymm13
               vfmadd231ps    C012, %ymm13, %ymm12
               vmovaps        %ymm14, 1056(%rdi)
               vmovaps        1056(%rsi), %ymm1
               vmulps         C010, %ymm1, %ymm0
               vmovaps        640(%rsi), %ymm3
               vmulps         C025, %ymm3, %ymm2
               vmovaps        384(%rsi), %ymm5
               vmulps         C012, %ymm5, %ymm4
               vmovaps        %ymm6, 1344(%rdi)
               vmovaps        512(%rsi), %ymm9
               vfmadd231ps    C020, %ymm9, %ymm8
               vmovaps        736(%rsi), %ymm11
               vmulps         C028, %ymm11, %ymm10
               vmovaps        1088(%rsi), %ymm13
               vfmadd231ps    C020, %ymm13, %ymm12
               vmovaps        416(%rsi), %ymm15
               vfmadd231ps    C027, %ymm15, %ymm14
               vmovaps        1120(%rsi), %ymm1
               vmulps         C024, %ymm1, %ymm0
               vmovaps        704(%rsi), %ymm3
               vmulps         C012, %ymm3, %ymm2
               vmovaps        352(%rsi), %ymm5
               vfmadd231ps    C026, %ymm5, %ymm4
               vmovaps        480(%rsi), %ymm7
               vfmadd231ps    C028, %ymm7, %ymm6
               vmovaps        %ymm8, 704(%rdi)
               vmovaps        %ymm10, 672(%rdi)
               vmovaps        %ymm12, 992(%rdi)
               vmovaps        %ymm0, 960(%rdi)
               vmovaps        768(%rsi), %ymm3
               vmulps         C020, %ymm3, %ymm2
               vmovaps        448(%rsi), %ymm5
               vmulps         C003, %ymm5, %ymm4
               vmovaps        %ymm6, 384(%rdi)
               vmovaps        288(%rsi), %ymm1
               vfmadd231ps    C028, %ymm1, %ymm0
               vmovaps        %ymm2, 640(%rdi)
               vmovaps        %ymm0, 416(%rdi)
               mov            Ldst, %rdi
               mov            Lwrk, %rsi
               vmovaps        (%rsi), %ymm0
               vmovaps        32(%rsi), %ymm1
               vmovaps        64(%rsi), %ymm2
               vmovaps        128(%rsi), %ymm3
               vmovaps        %ymm0, (%rdi)
               vmovaps        %ymm1, 32(%rdi)
               vmovaps        %ymm2, 64(%rdi)
               vmovaps        %ymm3, 128(%rdi)
               vmovups        R, %ymm0
               vmovaps        %ymm0, %ymm15
               vmulps         %ymm0, %ymm0, %ymm14
               vmulps         %ymm14, %ymm15, %ymm13
               vmulps         %ymm14, %ymm14, %ymm12
               vmulps         %ymm13, %ymm14, %ymm11
               vmulps         %ymm13, %ymm13, %ymm10
               vmulps         %ymm12, %ymm13, %ymm9
               vmulps         C003, %ymm14, %ymm14
               vmulps         C004, %ymm13, %ymm13
               vmulps         C005, %ymm12, %ymm12
               vmulps         C006, %ymm11, %ymm11
               vmulps         C007, %ymm10, %ymm10
               vmulps         C008, %ymm9, %ymm9
               vmovaps        0(%rsi), %ymm0
               vmovaps        64(%rsi), %ymm1
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vmovaps        192(%rsi), %ymm1
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vmovaps        384(%rsi), %ymm1
               vfmadd231ps    %ymm13, %ymm1, %ymm0
               vmovaps        640(%rsi), %ymm1
               vfmadd231ps    %ymm12, %ymm1, %ymm0
               vmovaps        960(%rsi), %ymm1
               vfmadd231ps    %ymm11, %ymm1, %ymm0
               vmovaps        1344(%rsi), %ymm1
               vfmadd231ps    %ymm10, %ymm1, %ymm0
               vmovaps        1792(%rsi), %ymm1
               vfmadd231ps    %ymm9, %ymm1, %ymm0
               vmovaps        %ymm0, 0(%rdi)
               vmovaps        96(%rsi), %ymm0
               vmovaps        64(%rsi), %ymm2
               vmovaps        32(%rsi), %ymm4
               vmovaps        224(%rsi), %ymm1
               vmovaps        192(%rsi), %ymm3
               vmovaps        160(%rsi), %ymm5
               vfmadd231ps    %ymm15, %ymm1, %ymm0
               vfmadd231ps    %ymm15, %ymm3, %ymm2
               vfmadd231ps    %ymm15, %ymm5, %ymm4
               vmovaps        416(%rsi), %ymm1
               vmovaps        384(%rsi), %ymm3
               vmovaps        352(%rsi), %ymm5
               vfmadd231ps    %ymm14, %ymm1, %ymm0
               vfmadd231ps    %ymm14, %ymm3, %ymm2
               vfmadd231ps    %ymm14, %ymm5, %ymm4
               vmovaps        672(%rsi), %ymm1
               vmovaps        640(%rsi), %ymm3
               vmovaps        608(%rsi), %ymm5
               vfmadd231ps    %ymm13, %ymm1, %ymm0
               vfmadd231ps    %ymm13, %ymm3, %ymm2
               vfmadd231ps    %ymm13, %ymm5, %ymm4
               vmovaps        992(%rsi), %ymm1
               vmovaps        960(%rsi), %ymm3
               vmovaps        928(%rsi), %ymm5
               vfmadd231ps    %ymm12, %ymm1, %ymm0
               vfmadd231ps    %ymm12, %ymm3, %ymm2
               vfmadd231ps    %ymm12, %ymm5, %ymm4
               vmovaps        1376(%rsi), %ymm1
               vmovaps        1344(%rsi), %ymm3
               vmovaps        1312(%rsi), %ymm5
               vfmadd231ps    %ymm11, %ymm1, %ymm0
               vfmadd231ps    %ymm11, %ymm3, %ymm2
               vfmadd231ps    %ymm11, %ymm5, %ymm4
               vmovaps        1824(%rsi), %ymm1
               vmovaps        1792(%rsi), %ymm3
               vmovaps        1760(%rsi), %ymm5
               vfmadd231ps    %ymm10, %ymm1, %ymm0
               vfmadd231ps    %ymm10, %ymm3, %ymm2
               vfmadd231ps    %ymm10, %ymm5, %ymm4
               vmovaps        %ymm0, 96(%rdi)
               vmovaps        %ymm2, 64(%rdi)
               vmovaps        %ymm4, 32(%rdi)
               vmovups        cosphi, %ymm0
               vmovups        sinphi, %ymm1
               mov            Lwrk, %rdi
               mov            Ldst, %rsi
               vaddps         64(%rsi), %ymm0, %ymm0
               vmovaps        32(%rsi), %ymm2
               vmovaps        96(%rsi), %ymm4
               vmovaps        %ymm0, 96(%rdi)
               vmovaps        %ymm2, 32(%rdi)
               vmovaps        %ymm4, 64(%rdi)
               mov            Lwrk, %rsi
               mov            Ldst, %rdi
               vmovaps        0(%r8), %ymm14
               vmovaps        0(%r9), %ymm15
               vmovaps        96(%rsi), %ymm2
               vmovaps        32(%rsi), %ymm3
               vmulps         %ymm3, %ymm14, %ymm0
               vmulps         %ymm2, %ymm15, %ymm1
               vfmadd231ps    %ymm2, %ymm14, %ymm0
               vfnmadd231ps   %ymm3, %ymm15, %ymm1
               vmovaps        %ymm0, 96(%rdi)
               vmovaps        %ymm1, 32(%rdi)
               mov            %rbp, %rsp
               pop            %rbp
               ret

               .align         32
C000:          .float         1.0000000000e-06  
               .float         1.0000000000e-06  
               .float         1.0000000000e-06  
               .float         1.0000000000e-06  
               .float         1.0000000000e-06  
               .float         1.0000000000e-06  
               .float         1.0000000000e-06  
               .float         1.0000000000e-06  
C001:          .float         1.0000000000e+00  
               .float         1.0000000000e+00  
               .float         1.0000000000e+00  
               .float         1.0000000000e+00  
               .float         1.0000000000e+00  
               .float         1.0000000000e+00  
               .float         1.0000000000e+00  
               .float         1.0000000000e+00  
C002:          .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
C003:          .float         5.0000000000e-01  
               .float         5.0000000000e-01  
               .float         5.0000000000e-01  
               .float         5.0000000000e-01  
               .float         5.0000000000e-01  
               .float         5.0000000000e-01  
               .float         5.0000000000e-01  
               .float         5.0000000000e-01  
C004:          .float         1.6666666667e-01  
               .float         1.6666666667e-01  
               .float         1.6666666667e-01  
               .float         1.6666666667e-01  
               .float         1.6666666667e-01  
               .float         1.6666666667e-01  
               .float         1.6666666667e-01  
               .float         1.6666666667e-01  
C005:          .float         4.1666666667e-02  
               .float         4.1666666667e-02  
               .float         4.1666666667e-02  
               .float         4.1666666667e-02  
               .float         4.1666666667e-02  
               .float         4.1666666667e-02  
               .float         4.1666666667e-02  
               .float         4.1666666667e-02  
C006:          .float         8.3333333333e-03  
               .float         8.3333333333e-03  
               .float         8.3333333333e-03  
               .float         8.3333333333e-03  
               .float         8.3333333333e-03  
               .float         8.3333333333e-03  
               .float         8.3333333333e-03  
               .float         8.3333333333e-03  
C007:          .float         1.3888888889e-03  
               .float         1.3888888889e-03  
               .float         1.3888888889e-03  
               .float         1.3888888889e-03  
               .float         1.3888888889e-03  
               .float         1.3888888889e-03  
               .float         1.3888888889e-03  
               .float         1.3888888889e-03  
C008:          .float         1.9841269841e-04  
               .float         1.9841269841e-04  
               .float         1.9841269841e-04  
               .float         1.9841269841e-04  
               .float         1.9841269841e-04  
               .float         1.9841269841e-04  
               .float         1.9841269841e-04  
               .float         1.9841269841e-04  
C009:          .float         -7.8125000000e-02 
               .float         -7.8125000000e-02 
               .float         -7.8125000000e-02 
               .float         -7.8125000000e-02 
               .float         -7.8125000000e-02 
               .float         -7.8125000000e-02 
               .float         -7.8125000000e-02 
               .float         -7.8125000000e-02 
C010:          .float         -3.1250000000e-01 
               .float         -3.1250000000e-01 
               .float         -3.1250000000e-01 
               .float         -3.1250000000e-01 
               .float         -3.1250000000e-01 
               .float         -3.1250000000e-01 
               .float         -3.1250000000e-01 
               .float         -3.1250000000e-01 
C011:          .float         -5.4687500000e-01 
               .float         -5.4687500000e-01 
               .float         -5.4687500000e-01 
               .float         -5.4687500000e-01 
               .float         -5.4687500000e-01 
               .float         -5.4687500000e-01 
               .float         -5.4687500000e-01 
               .float         -5.4687500000e-01 
C012:          .float         -5.0000000000e-01 
               .float         -5.0000000000e-01 
               .float         -5.0000000000e-01 
               .float         -5.0000000000e-01 
               .float         -5.0000000000e-01 
               .float         -5.0000000000e-01 
               .float         -5.0000000000e-01 
               .float         -5.0000000000e-01 
C013:          .float         1.4062500000e-01  
               .float         1.4062500000e-01  
               .float         1.4062500000e-01  
               .float         1.4062500000e-01  
               .float         1.4062500000e-01  
               .float         1.4062500000e-01  
               .float         1.4062500000e-01  
               .float         1.4062500000e-01  
C014:          .float         4.6875000000e-01  
               .float         4.6875000000e-01  
               .float         4.6875000000e-01  
               .float         4.6875000000e-01  
               .float         4.6875000000e-01  
               .float         4.6875000000e-01  
               .float         4.6875000000e-01  
               .float         4.6875000000e-01  
C015:          .float         3.2812500000e-01  
               .float         3.2812500000e-01  
               .float         3.2812500000e-01  
               .float         3.2812500000e-01  
               .float         3.2812500000e-01  
               .float         3.2812500000e-01  
               .float         3.2812500000e-01  
               .float         3.2812500000e-01  
C016:          .float         1.5625000000e-01  
               .float         1.5625000000e-01  
               .float         1.5625000000e-01  
               .float         1.5625000000e-01  
               .float         1.5625000000e-01  
               .float         1.5625000000e-01  
               .float         1.5625000000e-01  
               .float         1.5625000000e-01  
C017:          .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
C018:          .float         -1.8750000000e-01 
               .float         -1.8750000000e-01 
               .float         -1.8750000000e-01 
               .float         -1.8750000000e-01 
               .float         -1.8750000000e-01 
               .float         -1.8750000000e-01 
               .float         -1.8750000000e-01 
               .float         -1.8750000000e-01 
C019:          .float         -1.0937500000e-01 
               .float         -1.0937500000e-01 
               .float         -1.0937500000e-01 
               .float         -1.0937500000e-01 
               .float         -1.0937500000e-01 
               .float         -1.0937500000e-01 
               .float         -1.0937500000e-01 
               .float         -1.0937500000e-01 
C020:          .float         1.2500000000e-01  
               .float         1.2500000000e-01  
               .float         1.2500000000e-01  
               .float         1.2500000000e-01  
               .float         1.2500000000e-01  
               .float         1.2500000000e-01  
               .float         1.2500000000e-01  
               .float         1.2500000000e-01  
C021:          .float         1.5625000000e-02  
               .float         1.5625000000e-02  
               .float         1.5625000000e-02  
               .float         1.5625000000e-02  
               .float         1.5625000000e-02  
               .float         1.5625000000e-02  
               .float         1.5625000000e-02  
               .float         1.5625000000e-02  
C022:          .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
C023:          .float         3.1250000000e-02  
               .float         3.1250000000e-02  
               .float         3.1250000000e-02  
               .float         3.1250000000e-02  
               .float         3.1250000000e-02  
               .float         3.1250000000e-02  
               .float         3.1250000000e-02  
               .float         3.1250000000e-02  
C024:          .float         6.2500000000e-02  
               .float         6.2500000000e-02  
               .float         6.2500000000e-02  
               .float         6.2500000000e-02  
               .float         6.2500000000e-02  
               .float         6.2500000000e-02  
               .float         6.2500000000e-02  
               .float         6.2500000000e-02  
C025:          .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
C026:          .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
C027:          .float         -7.5000000000e-01 
               .float         -7.5000000000e-01 
               .float         -7.5000000000e-01 
               .float         -7.5000000000e-01 
               .float         -7.5000000000e-01 
               .float         -7.5000000000e-01 
               .float         -7.5000000000e-01 
               .float         -7.5000000000e-01 
C028:          .float         2.5000000000e-01  
               .float         2.5000000000e-01  
               .float         2.5000000000e-01  
               .float         2.5000000000e-01  
               .float         2.5000000000e-01  
               .float         2.5000000000e-01  
               .float         2.5000000000e-01  
               .float         2.5000000000e-01  
ABS:           .long          0x7fffffff
               .long          0x7fffffff
               .long          0x7fffffff
               .long          0x7fffffff
               .long          0x7fffffff
               .long          0x7fffffff
               .long          0x7fffffff
               .long          0x7fffffff
