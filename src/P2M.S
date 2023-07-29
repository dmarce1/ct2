#define        x              %ymm0
#define        y              %ymm1
#define        z              %ymm2
#define        r2             %ymm3
#define        ax0            %ymm4
#define        ay0            %ymm5
#define        M              %rdi
#define        Y              %rsi

               .global        M2M

               .text

M2M:           push           %rbp
               mov            %rsp, %rbp
               sub            $1568, %rsp
               and            $0xffffffffffffffc0, %rsp
               mov            %rsp, %rsi
               vmovaps        C000, %ymm9
               vmovaps        %ymm9, (Y)
               vmovaps        %ymm9, %ymm8
               vmulps         C001, x, %ymm8
               vmulps         C001, x, %ymm9
               vmovaps        %ymm8, 96(Y)
               vmovaps        %ymm9, 32(Y)
               vmulps         C002, %ymm8, ax0
               vmulps         C002, %ymm9, ay0
               vmulps         y, ay0, %ymm8
               vmulps         x, ay0, %ymm9
               vfmsub231ps    x, ax0, %ymm8
               vfmadd231ps    y, ax0, %ymm9
               vmovaps        %ymm8, 256(Y)
               vmovaps        %ymm9, 128(Y)
               vmulps         C003, %ymm8, ax0
               vmulps         C003, %ymm9, ay0
               vmulps         y, ay0, %ymm8
               vmulps         x, ay0, %ymm9
               vfmsub231ps    x, ax0, %ymm8
               vfmadd231ps    y, ax0, %ymm9
               vmovaps        %ymm8, 480(Y)
               vmovaps        %ymm9, 288(Y)
               vmulps         C004, %ymm8, ax0
               vmulps         C004, %ymm9, ay0
               vmulps         y, ay0, %ymm8
               vmulps         x, ay0, %ymm9
               vfmsub231ps    x, ax0, %ymm8
               vfmadd231ps    y, ax0, %ymm9
               vmovaps        %ymm8, 768(Y)
               vmovaps        %ymm9, 512(Y)
               vmulps         C005, %ymm8, ax0
               vmulps         C005, %ymm9, ay0
               vmulps         y, ay0, %ymm8
               vmulps         x, ay0, %ymm9
               vfmsub231ps    x, ax0, %ymm8
               vfmadd231ps    y, ax0, %ymm9
               vmovaps        %ymm8, 1120(Y)
               vmovaps        %ymm9, 800(Y)
               vmulps         C006, %ymm8, ax0
               vmulps         C006, %ymm9, ay0
               vmulps         y, ay0, %ymm8
               vmulps         x, ay0, %ymm9
               vfmsub231ps    x, ax0, %ymm8
               vfmadd231ps    y, ax0, %ymm9
               vmovaps        %ymm8, 1536(Y)
               vmovaps        %ymm9, 1152(Y)
               vmulps         C007, r2, %ymm9
               vmovaps        %ymm9, 192(Y)
               vmovaps        z, 64(Y)
               vmulps         C008, r2, %ymm9
               vmulps         64(Y), %ymm9, %ymm9
               vmovaps        %ymm9, 384(Y)
               vmulps         C009, r2, ax0
               vmulps         96(Y), ax0, %ymm8
               vmulps         1(Y), ax0, %ymm9
               vmovaps        %ymm8, 416(Y)
               vmovaps        %ymm9, 352(Y)
               vmulps         96(Y), z, %ymm8
               vmulps         32(Y), z, %ymm9
               vmovaps        %ymm8, 416(Y)
               vmovaps        %ymm9, 352(Y)
               vmulps         C010, r2, %ymm9
               vmulps         192(Y), %ymm9, %ymm9
               vmovaps        %ymm9, 640(Y)
               vmulps         C011, r2, ax0
               vmulps         224(Y), ax0, %ymm8
               vmulps         5(Y), ax0, %ymm9
               vmovaps        %ymm8, 672(Y)
               vmovaps        %ymm9, 608(Y)
               vmulps         C012, z, ay0
               vmovaps        416(Y), %ymm8
               vmovaps        352(Y), %ymm9
               vfmadd231ps    224(Y), ay0, %ymm8
               vfmadd231ps    160(Y), ay0, %ymm9
               vmovaps        %ymm8, 416(Y)
               vmovaps        %ymm9, 352(Y)
               vmulps         C013, r2, ax0
               vmulps         256(Y), ax0, %ymm8
               vmulps         4(Y), ax0, %ymm9
               vmovaps        %ymm8, 704(Y)
               vmovaps        %ymm9, 576(Y)
               vmulps         256(Y), z, %ymm8
               vmulps         128(Y), z, %ymm9
               vmovaps        %ymm8, 704(Y)
               vmovaps        %ymm9, 576(Y)
               vmulps         C014, r2, %ymm9
               vmulps         384(Y), %ymm9, %ymm9
               vmovaps        %ymm9, 960(Y)
               vmulps         C015, r2, ax0
               vmulps         416(Y), ax0, %ymm8
               vmulps         11(Y), ax0, %ymm9
               vmovaps        %ymm8, 992(Y)
               vmovaps        %ymm9, 928(Y)
               vmulps         C016, z, ay0
               vmovaps        672(Y), %ymm8
               vmovaps        608(Y), %ymm9
               vfmadd231ps    416(Y), ay0, %ymm8
               vfmadd231ps    352(Y), ay0, %ymm9
               vmovaps        %ymm8, 672(Y)
               vmovaps        %ymm9, 608(Y)
               vmulps         C017, r2, ax0
               vmulps         448(Y), ax0, %ymm8
               vmulps         10(Y), ax0, %ymm9
               vmovaps        %ymm8, 1024(Y)
               vmovaps        %ymm9, 896(Y)
               vmulps         C018, z, ay0
               vmovaps        704(Y), %ymm8
               vmovaps        576(Y), %ymm9
               vfmadd231ps    448(Y), ay0, %ymm8
               vfmadd231ps    320(Y), ay0, %ymm9
               vmovaps        %ymm8, 704(Y)
               vmovaps        %ymm9, 576(Y)
               vmulps         C010, r2, ax0
               vmulps         480(Y), ax0, %ymm8
               vmulps         9(Y), ax0, %ymm9
               vmovaps        %ymm8, 1056(Y)
               vmovaps        %ymm9, 864(Y)
               vmulps         480(Y), z, %ymm8
               vmulps         288(Y), z, %ymm9
               vmovaps        %ymm8, 1056(Y)
               vmovaps        %ymm9, 864(Y)
               vmulps         C019, r2, %ymm9
               vmulps         640(Y), %ymm9, %ymm9
               vmovaps        %ymm9, 1344(Y)
               vmulps         C020, r2, ax0
               vmulps         672(Y), ax0, %ymm8
               vmulps         19(Y), ax0, %ymm9
               vmovaps        %ymm8, 1376(Y)
               vmovaps        %ymm9, 1312(Y)
               vmulps         C021, z, ay0
               vmovaps        992(Y), %ymm8
               vmovaps        928(Y), %ymm9
               vfmadd231ps    672(Y), ay0, %ymm8
               vfmadd231ps    608(Y), ay0, %ymm9
               vmovaps        %ymm8, 992(Y)
               vmovaps        %ymm9, 928(Y)
               vmulps         C022, r2, ax0
               vmulps         704(Y), ax0, %ymm8
               vmulps         18(Y), ax0, %ymm9
               vmovaps        %ymm8, 1408(Y)
               vmovaps        %ymm9, 1280(Y)
               vmulps         C023, z, ay0
               vmovaps        1024(Y), %ymm8
               vmovaps        896(Y), %ymm9
               vfmadd231ps    704(Y), ay0, %ymm8
               vfmadd231ps    576(Y), ay0, %ymm9
               vmovaps        %ymm8, 1024(Y)
               vmovaps        %ymm9, 896(Y)
               vmulps         C024, r2, ax0
               vmulps         736(Y), ax0, %ymm8
               vmulps         17(Y), ax0, %ymm9
               vmovaps        %ymm8, 1440(Y)
               vmovaps        %ymm9, 1248(Y)
               vmulps         C025, z, ay0
               vmovaps        1056(Y), %ymm8
               vmovaps        864(Y), %ymm9
               vfmadd231ps    736(Y), ay0, %ymm8
               vfmadd231ps    544(Y), ay0, %ymm9
               vmovaps        %ymm8, 1056(Y)
               vmovaps        %ymm9, 864(Y)
               vmulps         C026, r2, ax0
               vmulps         768(Y), ax0, %ymm8
               vmulps         16(Y), ax0, %ymm9
               vmovaps        %ymm8, 1472(Y)
               vmovaps        %ymm9, 1216(Y)
               vmulps         768(Y), z, %ymm8
               vmulps         512(Y), z, %ymm9
               vmovaps        %ymm8, 1472(Y)
               vmovaps        %ymm9, 1216(Y)
               vmulps         C027, z, ay0
               vmovaps        1376(Y), %ymm8
               vmovaps        1312(Y), %ymm9
               vfmadd231ps    992(Y), ay0, %ymm8
               vfmadd231ps    928(Y), ay0, %ymm9
               vmovaps        %ymm8, 1376(Y)
               vmovaps        %ymm9, 1312(Y)
               vmulps         C028, z, ay0
               vmovaps        1408(Y), %ymm8
               vmovaps        1280(Y), %ymm9
               vfmadd231ps    1024(Y), ay0, %ymm8
               vfmadd231ps    896(Y), ay0, %ymm9
               vmovaps        %ymm8, 1408(Y)
               vmovaps        %ymm9, 1280(Y)
               vmulps         C029, z, ay0
               vmovaps        1440(Y), %ymm8
               vmovaps        1248(Y), %ymm9
               vfmadd231ps    1056(Y), ay0, %ymm8
               vfmadd231ps    864(Y), ay0, %ymm9
               vmovaps        %ymm8, 1440(Y)
               vmovaps        %ymm9, 1248(Y)
               vmulps         C030, z, ay0
               vmovaps        1472(Y), %ymm8
               vmovaps        1216(Y), %ymm9
               vfmadd231ps    1088(Y), ay0, %ymm8
               vfmadd231ps    832(Y), ay0, %ymm9
               vmovaps        %ymm8, 1472(Y)
               vmovaps        %ymm9, 1216(Y)
               vmulps         1120(Y), z, %ymm8
               vmulps         800(Y), z, %ymm9
               vmovaps        %ymm8, 1952(Y)
               vmovaps        %ymm9, 1632(Y)
               vmovaps        0(M), %ymm0
               vmovaps        32(M), %ymm1
               vmovaps        64(M), %ymm2
               vmovaps        96(M), %ymm3
               vmovaps        128(M), %ymm4
               vmovaps        160(M), %ymm5
               vmovaps        192(M), %ymm6
               vmovaps        224(M), %ymm7
               vmovaps        256(M), %ymm8
               vmovaps        288(M), %ymm9
               vmovaps        320(M), %ymm10
               vmovaps        352(M), %ymm11
               vmovaps        384(M), %ymm12
               vmovaps        416(M), %ymm13
               vmovaps        448(M), %ymm14
               vmovaps        480(M), %ymm15
               vaddps         0(Y), %ymm0, %ymm0
               vaddps         32(Y), %ymm1, %ymm1
               vaddps         64(Y), %ymm2, %ymm2
               vaddps         96(Y), %ymm3, %ymm3
               vaddps         128(Y), %ymm4, %ymm4
               vaddps         160(Y), %ymm5, %ymm5
               vaddps         192(Y), %ymm6, %ymm6
               vaddps         224(Y), %ymm7, %ymm7
               vaddps         256(Y), %ymm8, %ymm8
               vaddps         288(Y), %ymm9, %ymm9
               vaddps         320(Y), %ymm10, %ymm10
               vaddps         352(Y), %ymm11, %ymm11
               vaddps         384(Y), %ymm12, %ymm12
               vaddps         416(Y), %ymm13, %ymm13
               vaddps         448(Y), %ymm14, %ymm14
               vaddps         480(Y), %ymm15, %ymm15
               vmovaps        %ymm0, 0(M)
               vmovaps        %ymm1, 32(M)
               vmovaps        %ymm2, 64(M)
               vmovaps        %ymm3, 96(M)
               vmovaps        %ymm4, 128(M)
               vmovaps        %ymm5, 160(M)
               vmovaps        %ymm6, 192(M)
               vmovaps        %ymm7, 224(M)
               vmovaps        %ymm8, 256(M)
               vmovaps        %ymm9, 288(M)
               vmovaps        %ymm10, 320(M)
               vmovaps        %ymm11, 352(M)
               vmovaps        %ymm12, 384(M)
               vmovaps        %ymm13, 416(M)
               vmovaps        %ymm14, 448(M)
               vmovaps        %ymm15, 480(M)
               vmovaps        512(M), %ymm0
               vmovaps        544(M), %ymm1
               vmovaps        576(M), %ymm2
               vmovaps        608(M), %ymm3
               vmovaps        640(M), %ymm4
               vmovaps        672(M), %ymm5
               vmovaps        704(M), %ymm6
               vmovaps        736(M), %ymm7
               vmovaps        768(M), %ymm8
               vmovaps        800(M), %ymm9
               vmovaps        832(M), %ymm10
               vmovaps        864(M), %ymm11
               vmovaps        896(M), %ymm12
               vmovaps        928(M), %ymm13
               vmovaps        960(M), %ymm14
               vmovaps        992(M), %ymm15
               vaddps         512(Y), %ymm0, %ymm0
               vaddps         544(Y), %ymm1, %ymm1
               vaddps         576(Y), %ymm2, %ymm2
               vaddps         608(Y), %ymm3, %ymm3
               vaddps         640(Y), %ymm4, %ymm4
               vaddps         672(Y), %ymm5, %ymm5
               vaddps         704(Y), %ymm6, %ymm6
               vaddps         736(Y), %ymm7, %ymm7
               vaddps         768(Y), %ymm8, %ymm8
               vaddps         800(Y), %ymm9, %ymm9
               vaddps         832(Y), %ymm10, %ymm10
               vaddps         864(Y), %ymm11, %ymm11
               vaddps         896(Y), %ymm12, %ymm12
               vaddps         928(Y), %ymm13, %ymm13
               vaddps         960(Y), %ymm14, %ymm14
               vaddps         992(Y), %ymm15, %ymm15
               vmovaps        %ymm0, 512(M)
               vmovaps        %ymm1, 544(M)
               vmovaps        %ymm2, 576(M)
               vmovaps        %ymm3, 608(M)
               vmovaps        %ymm4, 640(M)
               vmovaps        %ymm5, 672(M)
               vmovaps        %ymm6, 704(M)
               vmovaps        %ymm7, 736(M)
               vmovaps        %ymm8, 768(M)
               vmovaps        %ymm9, 800(M)
               vmovaps        %ymm10, 832(M)
               vmovaps        %ymm11, 864(M)
               vmovaps        %ymm12, 896(M)
               vmovaps        %ymm13, 928(M)
               vmovaps        %ymm14, 960(M)
               vmovaps        %ymm15, 992(M)
               vmovaps        1024(M), %ymm0
               vmovaps        1056(M), %ymm1
               vmovaps        1088(M), %ymm2
               vmovaps        1120(M), %ymm3
               vmovaps        1152(M), %ymm4
               vmovaps        1184(M), %ymm5
               vmovaps        1216(M), %ymm6
               vmovaps        1248(M), %ymm7
               vmovaps        1280(M), %ymm8
               vmovaps        1312(M), %ymm9
               vmovaps        1344(M), %ymm10
               vmovaps        1376(M), %ymm11
               vmovaps        1408(M), %ymm12
               vmovaps        1440(M), %ymm13
               vmovaps        1472(M), %ymm14
               vmovaps        1504(M), %ymm15
               vaddps         1024(Y), %ymm0, %ymm0
               vaddps         1056(Y), %ymm1, %ymm1
               vaddps         1088(Y), %ymm2, %ymm2
               vaddps         1120(Y), %ymm3, %ymm3
               vaddps         1152(Y), %ymm4, %ymm4
               vaddps         1184(Y), %ymm5, %ymm5
               vaddps         1216(Y), %ymm6, %ymm6
               vaddps         1248(Y), %ymm7, %ymm7
               vaddps         1280(Y), %ymm8, %ymm8
               vaddps         1312(Y), %ymm9, %ymm9
               vaddps         1344(Y), %ymm10, %ymm10
               vaddps         1376(Y), %ymm11, %ymm11
               vaddps         1408(Y), %ymm12, %ymm12
               vaddps         1440(Y), %ymm13, %ymm13
               vaddps         1472(Y), %ymm14, %ymm14
               vaddps         1504(Y), %ymm15, %ymm15
               vmovaps        %ymm0, 1024(M)
               vmovaps        %ymm1, 1056(M)
               vmovaps        %ymm2, 1088(M)
               vmovaps        %ymm3, 1120(M)
               vmovaps        %ymm4, 1152(M)
               vmovaps        %ymm5, 1184(M)
               vmovaps        %ymm6, 1216(M)
               vmovaps        %ymm7, 1248(M)
               vmovaps        %ymm8, 1280(M)
               vmovaps        %ymm9, 1312(M)
               vmovaps        %ymm10, 1344(M)
               vmovaps        %ymm11, 1376(M)
               vmovaps        %ymm12, 1408(M)
               vmovaps        %ymm13, 1440(M)
               vmovaps        %ymm14, 1472(M)
               vmovaps        %ymm15, 1504(M)
               vmovaps        1536(M), %ymm0
               vaddps         1536(Y), %ymm0, %ymm0
               vmovaps        %ymm0, 1536(M)
               mov            %rsp, %rbp
               pop            %rbp
               ret

               .align         32
C000:          .float         1.0000000000e+00  
               .float         1.0000000000e+00  
               .float         1.0000000000e+00  
               .float         1.0000000000e+00  
               .float         1.0000000000e+00  
               .float         1.0000000000e+00  
               .float         1.0000000000e+00  
               .float         1.0000000000e+00  
C001:          .float         5.0000000000e-01  
               .float         5.0000000000e-01  
               .float         5.0000000000e-01  
               .float         5.0000000000e-01  
               .float         5.0000000000e-01  
               .float         5.0000000000e-01  
               .float         5.0000000000e-01  
               .float         5.0000000000e-01  
C002:          .float         2.5000000000e-01  
               .float         2.5000000000e-01  
               .float         2.5000000000e-01  
               .float         2.5000000000e-01  
               .float         2.5000000000e-01  
               .float         2.5000000000e-01  
               .float         2.5000000000e-01  
               .float         2.5000000000e-01  
C003:          .float         1.6666666667e-01  
               .float         1.6666666667e-01  
               .float         1.6666666667e-01  
               .float         1.6666666667e-01  
               .float         1.6666666667e-01  
               .float         1.6666666667e-01  
               .float         1.6666666667e-01  
               .float         1.6666666667e-01  
C004:          .float         1.2500000000e-01  
               .float         1.2500000000e-01  
               .float         1.2500000000e-01  
               .float         1.2500000000e-01  
               .float         1.2500000000e-01  
               .float         1.2500000000e-01  
               .float         1.2500000000e-01  
               .float         1.2500000000e-01  
C005:          .float         1.0000000000e-01  
               .float         1.0000000000e-01  
               .float         1.0000000000e-01  
               .float         1.0000000000e-01  
               .float         1.0000000000e-01  
               .float         1.0000000000e-01  
               .float         1.0000000000e-01  
               .float         1.0000000000e-01  
C006:          .float         8.3333333333e-02  
               .float         8.3333333333e-02  
               .float         8.3333333333e-02  
               .float         8.3333333333e-02  
               .float         8.3333333333e-02  
               .float         8.3333333333e-02  
               .float         8.3333333333e-02  
               .float         8.3333333333e-02  
C007:          .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
C008:          .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
C009:          .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
C010:          .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
C011:          .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
C012:          .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
C013:          .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
C014:          .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
C015:          .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
C016:          .float         4.6666666667e-01  
               .float         4.6666666667e-01  
               .float         4.6666666667e-01  
               .float         4.6666666667e-01  
               .float         4.6666666667e-01  
               .float         4.6666666667e-01  
               .float         4.6666666667e-01  
               .float         4.6666666667e-01  
C017:          .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
C018:          .float         5.8333333333e-01  
               .float         5.8333333333e-01  
               .float         5.8333333333e-01  
               .float         5.8333333333e-01  
               .float         5.8333333333e-01  
               .float         5.8333333333e-01  
               .float         5.8333333333e-01  
               .float         5.8333333333e-01  
C019:          .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
C020:          .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
C021:          .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
C022:          .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
C023:          .float         4.2857142857e-01  
               .float         4.2857142857e-01  
               .float         4.2857142857e-01  
               .float         4.2857142857e-01  
               .float         4.2857142857e-01  
               .float         4.2857142857e-01  
               .float         4.2857142857e-01  
               .float         4.2857142857e-01  
C024:          .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
C025:          .float         5.6250000000e-01  
               .float         5.6250000000e-01  
               .float         5.6250000000e-01  
               .float         5.6250000000e-01  
               .float         5.6250000000e-01  
               .float         5.6250000000e-01  
               .float         5.6250000000e-01  
               .float         5.6250000000e-01  
C026:          .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
C027:          .float         3.1428571429e-01  
               .float         3.1428571429e-01  
               .float         3.1428571429e-01  
               .float         3.1428571429e-01  
               .float         3.1428571429e-01  
               .float         3.1428571429e-01  
               .float         3.1428571429e-01  
               .float         3.1428571429e-01  
C028:          .float         3.4375000000e-01  
               .float         3.4375000000e-01  
               .float         3.4375000000e-01  
               .float         3.4375000000e-01  
               .float         3.4375000000e-01  
               .float         3.4375000000e-01  
               .float         3.4375000000e-01  
               .float         3.4375000000e-01  
C029:          .float         4.0740740741e-01  
               .float         4.0740740741e-01  
               .float         4.0740740741e-01  
               .float         4.0740740741e-01  
               .float         4.0740740741e-01  
               .float         4.0740740741e-01  
               .float         4.0740740741e-01  
               .float         4.0740740741e-01  
C030:          .float         5.5000000000e-01  
               .float         5.5000000000e-01  
               .float         5.5000000000e-01  
               .float         5.5000000000e-01  
               .float         5.5000000000e-01  
               .float         5.5000000000e-01  
               .float         5.5000000000e-01  
               .float         5.5000000000e-01  
