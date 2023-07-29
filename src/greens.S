#define        x              %ymm0
#define        y              %ymm1
#define        z              %ymm2
#define        r2             %ymm3
#define        r2inv          %ymm4
#define        ax0            %ymm5
#define        ay0            %ymm6
#define        Optr           %rdi

               .global        greens

               .text

greens:        vmulps         z, z, r2
               vfmadd231ps    y, y, r2
               vfmadd231ps    x, x, r2
               vmovaps        C000, %ymm8
               vdivps         r2, %ymm8, r2inv
               vrsqrtps       r2, %ymm8
               vmulps         C000, %ymm8, %ymm8
               vmovaps        %ymm8, (Optr)
               vmulps         C000, r2inv, %ymm9
               vmulps         %ymm9, x, x
               vmulps         %ymm9, y, y
               vmulps         %ymm9, z, z
               vmulps         %ymm8, x, %ymm9
               vmulps         %ymm8, y, %ymm8
               vmovaps        %ymm9, 96(Optr)
               vmovaps        %ymm8, 32(Optr)
               vmulps         C001, %ymm9, ax0
               vmulps         C001, %ymm8, ay0
               vmulps         y, ay0, %ymm9
               vmulps         x, ay0, %ymm8
               vfmsub231ps    x, ax0, %ymm9
               vfmadd231ps    y, ax0, %ymm8
               vmovaps        %ymm9, 256(Optr)
               vmovaps        %ymm8, 128(Optr)
               vmulps         C002, %ymm9, ax0
               vmulps         C002, %ymm8, ay0
               vmulps         y, ay0, %ymm9
               vmulps         x, ay0, %ymm8
               vfmsub231ps    x, ax0, %ymm9
               vfmadd231ps    y, ax0, %ymm8
               vmovaps        %ymm9, 480(Optr)
               vmovaps        %ymm8, 288(Optr)
               vmulps         C003, %ymm9, ax0
               vmulps         C003, %ymm8, ay0
               vmulps         y, ay0, %ymm9
               vmulps         x, ay0, %ymm8
               vfmsub231ps    x, ax0, %ymm9
               vfmadd231ps    y, ax0, %ymm8
               vmovaps        %ymm9, 768(Optr)
               vmovaps        %ymm8, 512(Optr)
               vmulps         C004, %ymm9, ax0
               vmulps         C004, %ymm8, ay0
               vmulps         y, ay0, %ymm9
               vmulps         x, ay0, %ymm8
               vfmsub231ps    x, ax0, %ymm9
               vfmadd231ps    y, ax0, %ymm8
               vmovaps        %ymm9, 1120(Optr)
               vmovaps        %ymm8, 800(Optr)
               vmulps         C005, %ymm9, ax0
               vmulps         C005, %ymm8, ay0
               vmulps         y, ay0, %ymm9
               vmulps         x, ay0, %ymm8
               vfmsub231ps    x, ax0, %ymm9
               vfmadd231ps    y, ax0, %ymm8
               vmovaps        %ymm9, 1536(Optr)
               vmovaps        %ymm8, 1152(Optr)
               vmulps         C006, %ymm9, ax0
               vmulps         C006, %ymm8, ay0
               vmulps         y, ay0, %ymm9
               vmulps         x, ay0, %ymm8
               vfmsub231ps    x, ax0, %ymm9
               vfmadd231ps    y, ax0, %ymm8
               vmovaps        %ymm9, 2016(Optr)
               vmovaps        %ymm8, 1568(Optr)
               vmulps         (Optr), r2inv, %ymm8
               vmulps         C000, %ymm8, %ymm8
               vmovaps        %ymm8, 192(Optr)
               vmovaps        %ymm9, 64(Optr)
               vmulps         C001, z, ax0
               vmulps         32(Optr), ax0, %ymm8
               vfmadd231ps    64(Optr), ax0, %ymm8
               vfmadd231ps    96(Optr), ax0, %ymm8
               vmovaps        %ymm8, 160(Optr)
               vmulps         C007, r2inv, %ymm8
               vmulps         64(Optr), %ymm8, %ymm8
               vmovaps        %ymm8, 384(Optr)
               vmulps         C008, r2inv, ax0
               vmulps         96(Optr), ax0, %ymm9
               vmulps         32(Optr), ax0, %ymm8
               vmovaps        %ymm9, 416(Optr)
               vmovaps        %ymm8, 352(Optr)
               vmulps         C002, z, ax0
               vmulps         128(Optr), ax0, %ymm8
               vfmadd231ps    160(Optr), ax0, %ymm8
               vfmadd231ps    192(Optr), ax0, %ymm8
               vfmadd231ps    224(Optr), ax0, %ymm8
               vfmadd231ps    256(Optr), ax0, %ymm8
               vmovaps        %ymm8, 320(Optr)
               vmulps         C009, r2inv, %ymm8
               vmulps         192(Optr), %ymm8, %ymm8
               vmovaps        %ymm8, 640(Optr)
               vmulps         C010, r2inv, ax0
               vmulps         224(Optr), ax0, %ymm9
               vmulps         160(Optr), ax0, %ymm8
               vmovaps        %ymm9, 672(Optr)
               vmovaps        %ymm8, 608(Optr)
               vmulps         C011, r2inv, ax0
               vmulps         256(Optr), ax0, %ymm9
               vmulps         128(Optr), ax0, %ymm8
               vmovaps        %ymm9, 704(Optr)
               vmovaps        %ymm8, 576(Optr)
               vmulps         C003, z, ax0
               vmulps         288(Optr), ax0, %ymm8
               vfmadd231ps    320(Optr), ax0, %ymm8
               vfmadd231ps    352(Optr), ax0, %ymm8
               vfmadd231ps    384(Optr), ax0, %ymm8
               vfmadd231ps    416(Optr), ax0, %ymm8
               vfmadd231ps    448(Optr), ax0, %ymm8
               vfmadd231ps    480(Optr), ax0, %ymm8
               vmovaps        %ymm8, 544(Optr)
               vmulps         C012, r2inv, %ymm8
               vmulps         384(Optr), %ymm8, %ymm8
               vmovaps        %ymm8, 960(Optr)
               vmulps         C013, r2inv, ax0
               vmulps         416(Optr), ax0, %ymm9
               vmulps         352(Optr), ax0, %ymm8
               vmovaps        %ymm9, 992(Optr)
               vmovaps        %ymm8, 928(Optr)
               vmulps         C014, r2inv, ax0
               vmulps         448(Optr), ax0, %ymm9
               vmulps         320(Optr), ax0, %ymm8
               vmovaps        %ymm9, 1024(Optr)
               vmovaps        %ymm8, 896(Optr)
               vmulps         C015, r2inv, ax0
               vmulps         480(Optr), ax0, %ymm9
               vmulps         288(Optr), ax0, %ymm8
               vmovaps        %ymm9, 1056(Optr)
               vmovaps        %ymm8, 864(Optr)
               vmulps         C004, z, ax0
               vmulps         512(Optr), ax0, %ymm8
               vfmadd231ps    544(Optr), ax0, %ymm8
               vfmadd231ps    576(Optr), ax0, %ymm8
               vfmadd231ps    608(Optr), ax0, %ymm8
               vfmadd231ps    640(Optr), ax0, %ymm8
               vfmadd231ps    672(Optr), ax0, %ymm8
               vfmadd231ps    704(Optr), ax0, %ymm8
               vfmadd231ps    736(Optr), ax0, %ymm8
               vfmadd231ps    768(Optr), ax0, %ymm8
               vmovaps        %ymm8, 832(Optr)
               vmulps         C016, r2inv, %ymm8
               vmulps         640(Optr), %ymm8, %ymm8
               vmovaps        %ymm8, 1344(Optr)
               vmulps         C017, r2inv, ax0
               vmulps         672(Optr), ax0, %ymm9
               vmulps         608(Optr), ax0, %ymm8
               vmovaps        %ymm9, 1376(Optr)
               vmovaps        %ymm8, 1312(Optr)
               vmulps         C018, r2inv, ax0
               vmulps         704(Optr), ax0, %ymm9
               vmulps         576(Optr), ax0, %ymm8
               vmovaps        %ymm9, 1408(Optr)
               vmovaps        %ymm8, 1280(Optr)
               vmulps         C012, r2inv, ax0
               vmulps         736(Optr), ax0, %ymm9
               vmulps         544(Optr), ax0, %ymm8
               vmovaps        %ymm9, 1440(Optr)
               vmovaps        %ymm8, 1248(Optr)
               vmulps         C009, r2inv, ax0
               vmulps         768(Optr), ax0, %ymm9
               vmulps         512(Optr), ax0, %ymm8
               vmovaps        %ymm9, 1472(Optr)
               vmovaps        %ymm8, 1216(Optr)
               vmulps         C005, z, ax0
               vmulps         800(Optr), ax0, %ymm8
               vfmadd231ps    832(Optr), ax0, %ymm8
               vfmadd231ps    864(Optr), ax0, %ymm8
               vfmadd231ps    896(Optr), ax0, %ymm8
               vfmadd231ps    928(Optr), ax0, %ymm8
               vfmadd231ps    960(Optr), ax0, %ymm8
               vfmadd231ps    992(Optr), ax0, %ymm8
               vfmadd231ps    1024(Optr), ax0, %ymm8
               vfmadd231ps    1056(Optr), ax0, %ymm8
               vfmadd231ps    1088(Optr), ax0, %ymm8
               vfmadd231ps    1120(Optr), ax0, %ymm8
               vmovaps        %ymm8, 1184(Optr)
               vmulps         C019, r2inv, %ymm8
               vmulps         960(Optr), %ymm8, %ymm8
               vmovaps        %ymm8, 1792(Optr)
               vmulps         C020, r2inv, ax0
               vmulps         992(Optr), ax0, %ymm9
               vmulps         928(Optr), ax0, %ymm8
               vmovaps        %ymm9, 1824(Optr)
               vmovaps        %ymm8, 1760(Optr)
               vmulps         C021, r2inv, ax0
               vmulps         1024(Optr), ax0, %ymm9
               vmulps         896(Optr), ax0, %ymm8
               vmovaps        %ymm9, 1856(Optr)
               vmovaps        %ymm8, 1728(Optr)
               vmulps         C022, r2inv, ax0
               vmulps         1056(Optr), ax0, %ymm9
               vmulps         864(Optr), ax0, %ymm8
               vmovaps        %ymm9, 1888(Optr)
               vmovaps        %ymm8, 1696(Optr)
               vmulps         C023, r2inv, ax0
               vmulps         1088(Optr), ax0, %ymm9
               vmulps         832(Optr), ax0, %ymm8
               vmovaps        %ymm9, 1920(Optr)
               vmovaps        %ymm8, 1664(Optr)
               vmulps         C024, r2inv, ax0
               vmulps         1120(Optr), ax0, %ymm9
               vmulps         800(Optr), ax0, %ymm8
               vmovaps        %ymm9, 1952(Optr)
               vmovaps        %ymm8, 1632(Optr)
               vmulps         C006, z, ax0
               vmulps         1152(Optr), ax0, %ymm8
               vfmadd231ps    1184(Optr), ax0, %ymm8
               vfmadd231ps    1216(Optr), ax0, %ymm8
               vfmadd231ps    1248(Optr), ax0, %ymm8
               vfmadd231ps    1280(Optr), ax0, %ymm8
               vfmadd231ps    1312(Optr), ax0, %ymm8
               vfmadd231ps    1344(Optr), ax0, %ymm8
               vfmadd231ps    1376(Optr), ax0, %ymm8
               vfmadd231ps    1408(Optr), ax0, %ymm8
               vfmadd231ps    1440(Optr), ax0, %ymm8
               vfmadd231ps    1472(Optr), ax0, %ymm8
               vfmadd231ps    1504(Optr), ax0, %ymm8
               vfmadd231ps    1536(Optr), ax0, %ymm8
               vmovaps        %ymm8, 1600(Optr)
               ret

               .align         32
C000:          .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
               .float         -1.0000000000e+00 
C001:          .float         3.0000000000e+00  
               .float         3.0000000000e+00  
               .float         3.0000000000e+00  
               .float         3.0000000000e+00  
               .float         3.0000000000e+00  
               .float         3.0000000000e+00  
               .float         3.0000000000e+00  
               .float         3.0000000000e+00  
C002:          .float         5.0000000000e+00  
               .float         5.0000000000e+00  
               .float         5.0000000000e+00  
               .float         5.0000000000e+00  
               .float         5.0000000000e+00  
               .float         5.0000000000e+00  
               .float         5.0000000000e+00  
               .float         5.0000000000e+00  
C003:          .float         7.0000000000e+00  
               .float         7.0000000000e+00  
               .float         7.0000000000e+00  
               .float         7.0000000000e+00  
               .float         7.0000000000e+00  
               .float         7.0000000000e+00  
               .float         7.0000000000e+00  
               .float         7.0000000000e+00  
C004:          .float         9.0000000000e+00  
               .float         9.0000000000e+00  
               .float         9.0000000000e+00  
               .float         9.0000000000e+00  
               .float         9.0000000000e+00  
               .float         9.0000000000e+00  
               .float         9.0000000000e+00  
               .float         9.0000000000e+00  
C005:          .float         1.1000000000e+01  
               .float         1.1000000000e+01  
               .float         1.1000000000e+01  
               .float         1.1000000000e+01  
               .float         1.1000000000e+01  
               .float         1.1000000000e+01  
               .float         1.1000000000e+01  
               .float         1.1000000000e+01  
C006:          .float         1.3000000000e+01  
               .float         1.3000000000e+01  
               .float         1.3000000000e+01  
               .float         1.3000000000e+01  
               .float         1.3000000000e+01  
               .float         1.3000000000e+01  
               .float         1.3000000000e+01  
               .float         1.3000000000e+01  
C007:          .float         -4.0000000000e+00 
               .float         -4.0000000000e+00 
               .float         -4.0000000000e+00 
               .float         -4.0000000000e+00 
               .float         -4.0000000000e+00 
               .float         -4.0000000000e+00 
               .float         -4.0000000000e+00 
               .float         -4.0000000000e+00 
C008:          .float         -3.0000000000e+00 
               .float         -3.0000000000e+00 
               .float         -3.0000000000e+00 
               .float         -3.0000000000e+00 
               .float         -3.0000000000e+00 
               .float         -3.0000000000e+00 
               .float         -3.0000000000e+00 
               .float         -3.0000000000e+00 
C009:          .float         -9.0000000000e+00 
               .float         -9.0000000000e+00 
               .float         -9.0000000000e+00 
               .float         -9.0000000000e+00 
               .float         -9.0000000000e+00 
               .float         -9.0000000000e+00 
               .float         -9.0000000000e+00 
               .float         -9.0000000000e+00 
C010:          .float         -8.0000000000e+00 
               .float         -8.0000000000e+00 
               .float         -8.0000000000e+00 
               .float         -8.0000000000e+00 
               .float         -8.0000000000e+00 
               .float         -8.0000000000e+00 
               .float         -8.0000000000e+00 
               .float         -8.0000000000e+00 
C011:          .float         -5.0000000000e+00 
               .float         -5.0000000000e+00 
               .float         -5.0000000000e+00 
               .float         -5.0000000000e+00 
               .float         -5.0000000000e+00 
               .float         -5.0000000000e+00 
               .float         -5.0000000000e+00 
               .float         -5.0000000000e+00 
C012:          .float         -1.6000000000e+01 
               .float         -1.6000000000e+01 
               .float         -1.6000000000e+01 
               .float         -1.6000000000e+01 
               .float         -1.6000000000e+01 
               .float         -1.6000000000e+01 
               .float         -1.6000000000e+01 
               .float         -1.6000000000e+01 
C013:          .float         -1.5000000000e+01 
               .float         -1.5000000000e+01 
               .float         -1.5000000000e+01 
               .float         -1.5000000000e+01 
               .float         -1.5000000000e+01 
               .float         -1.5000000000e+01 
               .float         -1.5000000000e+01 
               .float         -1.5000000000e+01 
C014:          .float         -1.2000000000e+01 
               .float         -1.2000000000e+01 
               .float         -1.2000000000e+01 
               .float         -1.2000000000e+01 
               .float         -1.2000000000e+01 
               .float         -1.2000000000e+01 
               .float         -1.2000000000e+01 
               .float         -1.2000000000e+01 
C015:          .float         -7.0000000000e+00 
               .float         -7.0000000000e+00 
               .float         -7.0000000000e+00 
               .float         -7.0000000000e+00 
               .float         -7.0000000000e+00 
               .float         -7.0000000000e+00 
               .float         -7.0000000000e+00 
               .float         -7.0000000000e+00 
C016:          .float         -2.5000000000e+01 
               .float         -2.5000000000e+01 
               .float         -2.5000000000e+01 
               .float         -2.5000000000e+01 
               .float         -2.5000000000e+01 
               .float         -2.5000000000e+01 
               .float         -2.5000000000e+01 
               .float         -2.5000000000e+01 
C017:          .float         -2.4000000000e+01 
               .float         -2.4000000000e+01 
               .float         -2.4000000000e+01 
               .float         -2.4000000000e+01 
               .float         -2.4000000000e+01 
               .float         -2.4000000000e+01 
               .float         -2.4000000000e+01 
               .float         -2.4000000000e+01 
C018:          .float         -2.1000000000e+01 
               .float         -2.1000000000e+01 
               .float         -2.1000000000e+01 
               .float         -2.1000000000e+01 
               .float         -2.1000000000e+01 
               .float         -2.1000000000e+01 
               .float         -2.1000000000e+01 
               .float         -2.1000000000e+01 
C019:          .float         -3.6000000000e+01 
               .float         -3.6000000000e+01 
               .float         -3.6000000000e+01 
               .float         -3.6000000000e+01 
               .float         -3.6000000000e+01 
               .float         -3.6000000000e+01 
               .float         -3.6000000000e+01 
               .float         -3.6000000000e+01 
C020:          .float         -3.5000000000e+01 
               .float         -3.5000000000e+01 
               .float         -3.5000000000e+01 
               .float         -3.5000000000e+01 
               .float         -3.5000000000e+01 
               .float         -3.5000000000e+01 
               .float         -3.5000000000e+01 
               .float         -3.5000000000e+01 
C021:          .float         -3.2000000000e+01 
               .float         -3.2000000000e+01 
               .float         -3.2000000000e+01 
               .float         -3.2000000000e+01 
               .float         -3.2000000000e+01 
               .float         -3.2000000000e+01 
               .float         -3.2000000000e+01 
               .float         -3.2000000000e+01 
C022:          .float         -2.7000000000e+01 
               .float         -2.7000000000e+01 
               .float         -2.7000000000e+01 
               .float         -2.7000000000e+01 
               .float         -2.7000000000e+01 
               .float         -2.7000000000e+01 
               .float         -2.7000000000e+01 
               .float         -2.7000000000e+01 
C023:          .float         -2.0000000000e+01 
               .float         -2.0000000000e+01 
               .float         -2.0000000000e+01 
               .float         -2.0000000000e+01 
               .float         -2.0000000000e+01 
               .float         -2.0000000000e+01 
               .float         -2.0000000000e+01 
               .float         -2.0000000000e+01 
C024:          .float         -1.1000000000e+01 
               .float         -1.1000000000e+01 
               .float         -1.1000000000e+01 
               .float         -1.1000000000e+01 
               .float         -1.1000000000e+01 
               .float         -1.1000000000e+01 
               .float         -1.1000000000e+01 
               .float         -1.1000000000e+01 
