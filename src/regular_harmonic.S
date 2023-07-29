#define        x              %ymm0
#define        y              %ymm1
#define        z              %ymm2
#define        r2             %ymm3
#define        ax0            %ymm4
#define        ay0            %ymm5
#define        Y              %rdi

               .global        regular_harmonic

               .text

regular_harmonic:
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
               vmulps         C007, %ymm8, ax0
               vmulps         C007, %ymm9, ay0
               vmulps         y, ay0, %ymm8
               vmulps         x, ay0, %ymm9
               vfmsub231ps    x, ax0, %ymm8
               vfmadd231ps    y, ax0, %ymm9
               vmovaps        %ymm8, 2016(Y)
               vmovaps        %ymm9, 1568(Y)
               vmulps         C008, %ymm8, ax0
               vmulps         C008, %ymm9, ay0
               vmulps         y, ay0, %ymm8
               vmulps         x, ay0, %ymm9
               vfmsub231ps    x, ax0, %ymm8
               vfmadd231ps    y, ax0, %ymm9
               vmovaps        %ymm8, 2560(Y)
               vmovaps        %ymm9, 2048(Y)
               vmulps         C009, r2, %ymm9
               vmovaps        %ymm9, 192(Y)
               vmovaps        z, 64(Y)
               vmulps         C010, r2, %ymm9
               vmulps         64(Y), %ymm9, %ymm9
               vmovaps        %ymm9, 384(Y)
               vmulps         C011, r2, ax0
               vmulps         96(Y), ax0, %ymm8
               vmulps         1(Y), ax0, %ymm9
               vmovaps        %ymm8, 416(Y)
               vmovaps        %ymm9, 352(Y)
               vmulps         96(Y), z, %ymm8
               vmulps         32(Y), z, %ymm9
               vmovaps        %ymm8, 416(Y)
               vmovaps        %ymm9, 352(Y)
               vmulps         C012, r2, %ymm9
               vmulps         192(Y), %ymm9, %ymm9
               vmovaps        %ymm9, 640(Y)
               vmulps         C013, r2, ax0
               vmulps         224(Y), ax0, %ymm8
               vmulps         5(Y), ax0, %ymm9
               vmovaps        %ymm8, 672(Y)
               vmovaps        %ymm9, 608(Y)
               vmulps         C014, z, ay0
               vmovaps        416(Y), %ymm8
               vmovaps        352(Y), %ymm9
               vfmadd231ps    224(Y), ay0, %ymm8
               vfmadd231ps    160(Y), ay0, %ymm9
               vmovaps        %ymm8, 416(Y)
               vmovaps        %ymm9, 352(Y)
               vmulps         C015, r2, ax0
               vmulps         256(Y), ax0, %ymm8
               vmulps         4(Y), ax0, %ymm9
               vmovaps        %ymm8, 704(Y)
               vmovaps        %ymm9, 576(Y)
               vmulps         256(Y), z, %ymm8
               vmulps         128(Y), z, %ymm9
               vmovaps        %ymm8, 704(Y)
               vmovaps        %ymm9, 576(Y)
               vmulps         C016, r2, %ymm9
               vmulps         384(Y), %ymm9, %ymm9
               vmovaps        %ymm9, 960(Y)
               vmulps         C017, r2, ax0
               vmulps         416(Y), ax0, %ymm8
               vmulps         11(Y), ax0, %ymm9
               vmovaps        %ymm8, 992(Y)
               vmovaps        %ymm9, 928(Y)
               vmulps         C018, z, ay0
               vmovaps        672(Y), %ymm8
               vmovaps        608(Y), %ymm9
               vfmadd231ps    416(Y), ay0, %ymm8
               vfmadd231ps    352(Y), ay0, %ymm9
               vmovaps        %ymm8, 672(Y)
               vmovaps        %ymm9, 608(Y)
               vmulps         C019, r2, ax0
               vmulps         448(Y), ax0, %ymm8
               vmulps         10(Y), ax0, %ymm9
               vmovaps        %ymm8, 1024(Y)
               vmovaps        %ymm9, 896(Y)
               vmulps         C020, z, ay0
               vmovaps        704(Y), %ymm8
               vmovaps        576(Y), %ymm9
               vfmadd231ps    448(Y), ay0, %ymm8
               vfmadd231ps    320(Y), ay0, %ymm9
               vmovaps        %ymm8, 704(Y)
               vmovaps        %ymm9, 576(Y)
               vmulps         C012, r2, ax0
               vmulps         480(Y), ax0, %ymm8
               vmulps         9(Y), ax0, %ymm9
               vmovaps        %ymm8, 1056(Y)
               vmovaps        %ymm9, 864(Y)
               vmulps         480(Y), z, %ymm8
               vmulps         288(Y), z, %ymm9
               vmovaps        %ymm8, 1056(Y)
               vmovaps        %ymm9, 864(Y)
               vmulps         C021, r2, %ymm9
               vmulps         640(Y), %ymm9, %ymm9
               vmovaps        %ymm9, 1344(Y)
               vmulps         C022, r2, ax0
               vmulps         672(Y), ax0, %ymm8
               vmulps         19(Y), ax0, %ymm9
               vmovaps        %ymm8, 1376(Y)
               vmovaps        %ymm9, 1312(Y)
               vmulps         C023, z, ay0
               vmovaps        992(Y), %ymm8
               vmovaps        928(Y), %ymm9
               vfmadd231ps    672(Y), ay0, %ymm8
               vfmadd231ps    608(Y), ay0, %ymm9
               vmovaps        %ymm8, 992(Y)
               vmovaps        %ymm9, 928(Y)
               vmulps         C024, r2, ax0
               vmulps         704(Y), ax0, %ymm8
               vmulps         18(Y), ax0, %ymm9
               vmovaps        %ymm8, 1408(Y)
               vmovaps        %ymm9, 1280(Y)
               vmulps         C025, z, ay0
               vmovaps        1024(Y), %ymm8
               vmovaps        896(Y), %ymm9
               vfmadd231ps    704(Y), ay0, %ymm8
               vfmadd231ps    576(Y), ay0, %ymm9
               vmovaps        %ymm8, 1024(Y)
               vmovaps        %ymm9, 896(Y)
               vmulps         C026, r2, ax0
               vmulps         736(Y), ax0, %ymm8
               vmulps         17(Y), ax0, %ymm9
               vmovaps        %ymm8, 1440(Y)
               vmovaps        %ymm9, 1248(Y)
               vmulps         C027, z, ay0
               vmovaps        1056(Y), %ymm8
               vmovaps        864(Y), %ymm9
               vfmadd231ps    736(Y), ay0, %ymm8
               vfmadd231ps    544(Y), ay0, %ymm9
               vmovaps        %ymm8, 1056(Y)
               vmovaps        %ymm9, 864(Y)
               vmulps         C028, r2, ax0
               vmulps         768(Y), ax0, %ymm8
               vmulps         16(Y), ax0, %ymm9
               vmovaps        %ymm8, 1472(Y)
               vmovaps        %ymm9, 1216(Y)
               vmulps         768(Y), z, %ymm8
               vmulps         512(Y), z, %ymm9
               vmovaps        %ymm8, 1472(Y)
               vmovaps        %ymm9, 1216(Y)
               vmulps         C029, r2, %ymm9
               vmulps         960(Y), %ymm9, %ymm9
               vmovaps        %ymm9, 1792(Y)
               vmulps         C030, r2, ax0
               vmulps         992(Y), ax0, %ymm8
               vmulps         29(Y), ax0, %ymm9
               vmovaps        %ymm8, 1824(Y)
               vmovaps        %ymm9, 1760(Y)
               vmulps         C031, z, ay0
               vmovaps        1376(Y), %ymm8
               vmovaps        1312(Y), %ymm9
               vfmadd231ps    992(Y), ay0, %ymm8
               vfmadd231ps    928(Y), ay0, %ymm9
               vmovaps        %ymm8, 1376(Y)
               vmovaps        %ymm9, 1312(Y)
               vmulps         C032, r2, ax0
               vmulps         1024(Y), ax0, %ymm8
               vmulps         28(Y), ax0, %ymm9
               vmovaps        %ymm8, 1856(Y)
               vmovaps        %ymm9, 1728(Y)
               vmulps         C033, z, ay0
               vmovaps        1408(Y), %ymm8
               vmovaps        1280(Y), %ymm9
               vfmadd231ps    1024(Y), ay0, %ymm8
               vfmadd231ps    896(Y), ay0, %ymm9
               vmovaps        %ymm8, 1408(Y)
               vmovaps        %ymm9, 1280(Y)
               vmulps         C034, r2, ax0
               vmulps         1056(Y), ax0, %ymm8
               vmulps         27(Y), ax0, %ymm9
               vmovaps        %ymm8, 1888(Y)
               vmovaps        %ymm9, 1696(Y)
               vmulps         C035, z, ay0
               vmovaps        1440(Y), %ymm8
               vmovaps        1248(Y), %ymm9
               vfmadd231ps    1056(Y), ay0, %ymm8
               vfmadd231ps    864(Y), ay0, %ymm9
               vmovaps        %ymm8, 1440(Y)
               vmovaps        %ymm9, 1248(Y)
               vmulps         C036, r2, ax0
               vmulps         1088(Y), ax0, %ymm8
               vmulps         26(Y), ax0, %ymm9
               vmovaps        %ymm8, 1920(Y)
               vmovaps        %ymm9, 1664(Y)
               vmulps         C037, z, ay0
               vmovaps        1472(Y), %ymm8
               vmovaps        1216(Y), %ymm9
               vfmadd231ps    1088(Y), ay0, %ymm8
               vfmadd231ps    832(Y), ay0, %ymm9
               vmovaps        %ymm8, 1472(Y)
               vmovaps        %ymm9, 1216(Y)
               vmulps         C017, r2, ax0
               vmulps         1120(Y), ax0, %ymm8
               vmulps         25(Y), ax0, %ymm9
               vmovaps        %ymm8, 1952(Y)
               vmovaps        %ymm9, 1632(Y)
               vmulps         1120(Y), z, %ymm8
               vmulps         800(Y), z, %ymm9
               vmovaps        %ymm8, 1952(Y)
               vmovaps        %ymm9, 1632(Y)
               vmulps         C038, r2, %ymm9
               vmulps         1344(Y), %ymm9, %ymm9
               vmovaps        %ymm9, 2304(Y)
               vmulps         C039, r2, ax0
               vmulps         1376(Y), ax0, %ymm8
               vmulps         41(Y), ax0, %ymm9
               vmovaps        %ymm8, 2336(Y)
               vmovaps        %ymm9, 2272(Y)
               vmulps         C040, z, ay0
               vmovaps        1824(Y), %ymm8
               vmovaps        1760(Y), %ymm9
               vfmadd231ps    1376(Y), ay0, %ymm8
               vfmadd231ps    1312(Y), ay0, %ymm9
               vmovaps        %ymm8, 1824(Y)
               vmovaps        %ymm9, 1760(Y)
               vmulps         C041, r2, ax0
               vmulps         1408(Y), ax0, %ymm8
               vmulps         40(Y), ax0, %ymm9
               vmovaps        %ymm8, 2368(Y)
               vmovaps        %ymm9, 2240(Y)
               vmulps         C042, z, ay0
               vmovaps        1856(Y), %ymm8
               vmovaps        1728(Y), %ymm9
               vfmadd231ps    1408(Y), ay0, %ymm8
               vfmadd231ps    1280(Y), ay0, %ymm9
               vmovaps        %ymm8, 1856(Y)
               vmovaps        %ymm9, 1728(Y)
               vmulps         C043, r2, ax0
               vmulps         1440(Y), ax0, %ymm8
               vmulps         39(Y), ax0, %ymm9
               vmovaps        %ymm8, 2400(Y)
               vmovaps        %ymm9, 2208(Y)
               vmulps         C044, z, ay0
               vmovaps        1888(Y), %ymm8
               vmovaps        1696(Y), %ymm9
               vfmadd231ps    1440(Y), ay0, %ymm8
               vfmadd231ps    1248(Y), ay0, %ymm9
               vmovaps        %ymm8, 1888(Y)
               vmovaps        %ymm9, 1696(Y)
               vmulps         C030, r2, ax0
               vmulps         1472(Y), ax0, %ymm8
               vmulps         38(Y), ax0, %ymm9
               vmovaps        %ymm8, 2432(Y)
               vmovaps        %ymm9, 2176(Y)
               vmulps         C045, z, ay0
               vmovaps        1920(Y), %ymm8
               vmovaps        1664(Y), %ymm9
               vfmadd231ps    1472(Y), ay0, %ymm8
               vfmadd231ps    1216(Y), ay0, %ymm9
               vmovaps        %ymm8, 1920(Y)
               vmovaps        %ymm9, 1664(Y)
               vmulps         C046, r2, ax0
               vmulps         1504(Y), ax0, %ymm8
               vmulps         37(Y), ax0, %ymm9
               vmovaps        %ymm8, 2464(Y)
               vmovaps        %ymm9, 2144(Y)
               vmulps         C047, z, ay0
               vmovaps        1952(Y), %ymm8
               vmovaps        1632(Y), %ymm9
               vfmadd231ps    1504(Y), ay0, %ymm8
               vfmadd231ps    1184(Y), ay0, %ymm9
               vmovaps        %ymm8, 1952(Y)
               vmovaps        %ymm9, 1632(Y)
               vmulps         C048, r2, ax0
               vmulps         1536(Y), ax0, %ymm8
               vmulps         36(Y), ax0, %ymm9
               vmovaps        %ymm8, 2496(Y)
               vmovaps        %ymm9, 2112(Y)
               vmulps         1536(Y), z, %ymm8
               vmulps         1152(Y), z, %ymm9
               vmovaps        %ymm8, 2496(Y)
               vmovaps        %ymm9, 2112(Y)
               vmulps         C049, z, ay0
               vmovaps        2336(Y), %ymm8
               vmovaps        2272(Y), %ymm9
               vfmadd231ps    1824(Y), ay0, %ymm8
               vfmadd231ps    1760(Y), ay0, %ymm9
               vmovaps        %ymm8, 2336(Y)
               vmovaps        %ymm9, 2272(Y)
               vmulps         C002, z, ay0
               vmovaps        2368(Y), %ymm8
               vmovaps        2240(Y), %ymm9
               vfmadd231ps    1856(Y), ay0, %ymm8
               vfmadd231ps    1728(Y), ay0, %ymm9
               vmovaps        %ymm8, 2368(Y)
               vmovaps        %ymm9, 2240(Y)
               vmulps         C050, z, ay0
               vmovaps        2400(Y), %ymm8
               vmovaps        2208(Y), %ymm9
               vfmadd231ps    1888(Y), ay0, %ymm8
               vfmadd231ps    1696(Y), ay0, %ymm9
               vmovaps        %ymm8, 2400(Y)
               vmovaps        %ymm9, 2208(Y)
               vmulps         C051, z, ay0
               vmovaps        2432(Y), %ymm8
               vmovaps        2176(Y), %ymm9
               vfmadd231ps    1920(Y), ay0, %ymm8
               vfmadd231ps    1664(Y), ay0, %ymm9
               vmovaps        %ymm8, 2432(Y)
               vmovaps        %ymm9, 2176(Y)
               vmulps         C052, z, ay0
               vmovaps        2464(Y), %ymm8
               vmovaps        2144(Y), %ymm9
               vfmadd231ps    1952(Y), ay0, %ymm8
               vfmadd231ps    1632(Y), ay0, %ymm9
               vmovaps        %ymm8, 2464(Y)
               vmovaps        %ymm9, 2144(Y)
               vmulps         C053, z, ay0
               vmovaps        2496(Y), %ymm8
               vmovaps        2112(Y), %ymm9
               vfmadd231ps    1984(Y), ay0, %ymm8
               vfmadd231ps    1600(Y), ay0, %ymm9
               vmovaps        %ymm8, 2496(Y)
               vmovaps        %ymm9, 2112(Y)
               vmulps         2016(Y), z, %ymm8
               vmulps         1568(Y), z, %ymm9
               vmovaps        %ymm8, 3104(Y)
               vmovaps        %ymm9, 2656(Y)
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
C007:          .float         7.1428571429e-02  
               .float         7.1428571429e-02  
               .float         7.1428571429e-02  
               .float         7.1428571429e-02  
               .float         7.1428571429e-02  
               .float         7.1428571429e-02  
               .float         7.1428571429e-02  
               .float         7.1428571429e-02  
C008:          .float         6.2500000000e-02  
               .float         6.2500000000e-02  
               .float         6.2500000000e-02  
               .float         6.2500000000e-02  
               .float         6.2500000000e-02  
               .float         6.2500000000e-02  
               .float         6.2500000000e-02  
               .float         6.2500000000e-02  
C009:          .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
               .float         -2.5000000000e-01 
C010:          .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
               .float         -1.1111111111e-01 
C011:          .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
               .float         -1.2500000000e-01 
C012:          .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
               .float         -6.2500000000e-02 
C013:          .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
               .float         -6.6666666667e-02 
C014:          .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
               .float         6.2500000000e-01  
C015:          .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
               .float         -8.3333333333e-02 
C016:          .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
               .float         -4.0000000000e-02 
C017:          .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
               .float         -4.1666666667e-02 
C018:          .float         4.6666666667e-01  
               .float         4.6666666667e-01  
               .float         4.6666666667e-01  
               .float         4.6666666667e-01  
               .float         4.6666666667e-01  
               .float         4.6666666667e-01  
               .float         4.6666666667e-01  
               .float         4.6666666667e-01  
C019:          .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
               .float         -4.7619047619e-02 
C020:          .float         5.8333333333e-01  
               .float         5.8333333333e-01  
               .float         5.8333333333e-01  
               .float         5.8333333333e-01  
               .float         5.8333333333e-01  
               .float         5.8333333333e-01  
               .float         5.8333333333e-01  
               .float         5.8333333333e-01  
C021:          .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
               .float         -2.7777777778e-02 
C022:          .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
               .float         -2.8571428571e-02 
C023:          .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
               .float         3.7500000000e-01  
C024:          .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
               .float         -3.1250000000e-02 
C025:          .float         4.2857142857e-01  
               .float         4.2857142857e-01  
               .float         4.2857142857e-01  
               .float         4.2857142857e-01  
               .float         4.2857142857e-01  
               .float         4.2857142857e-01  
               .float         4.2857142857e-01  
               .float         4.2857142857e-01  
C026:          .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
               .float         -3.7037037037e-02 
C027:          .float         5.6250000000e-01  
               .float         5.6250000000e-01  
               .float         5.6250000000e-01  
               .float         5.6250000000e-01  
               .float         5.6250000000e-01  
               .float         5.6250000000e-01  
               .float         5.6250000000e-01  
               .float         5.6250000000e-01  
C028:          .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
               .float         -5.0000000000e-02 
C029:          .float         -2.0408163265e-02 
               .float         -2.0408163265e-02 
               .float         -2.0408163265e-02 
               .float         -2.0408163265e-02 
               .float         -2.0408163265e-02 
               .float         -2.0408163265e-02 
               .float         -2.0408163265e-02 
               .float         -2.0408163265e-02 
C030:          .float         -2.0833333333e-02 
               .float         -2.0833333333e-02 
               .float         -2.0833333333e-02 
               .float         -2.0833333333e-02 
               .float         -2.0833333333e-02 
               .float         -2.0833333333e-02 
               .float         -2.0833333333e-02 
               .float         -2.0833333333e-02 
C031:          .float         3.1428571429e-01  
               .float         3.1428571429e-01  
               .float         3.1428571429e-01  
               .float         3.1428571429e-01  
               .float         3.1428571429e-01  
               .float         3.1428571429e-01  
               .float         3.1428571429e-01  
               .float         3.1428571429e-01  
C032:          .float         -2.2222222222e-02 
               .float         -2.2222222222e-02 
               .float         -2.2222222222e-02 
               .float         -2.2222222222e-02 
               .float         -2.2222222222e-02 
               .float         -2.2222222222e-02 
               .float         -2.2222222222e-02 
               .float         -2.2222222222e-02 
C033:          .float         3.4375000000e-01  
               .float         3.4375000000e-01  
               .float         3.4375000000e-01  
               .float         3.4375000000e-01  
               .float         3.4375000000e-01  
               .float         3.4375000000e-01  
               .float         3.4375000000e-01  
               .float         3.4375000000e-01  
C034:          .float         -2.5000000000e-02 
               .float         -2.5000000000e-02 
               .float         -2.5000000000e-02 
               .float         -2.5000000000e-02 
               .float         -2.5000000000e-02 
               .float         -2.5000000000e-02 
               .float         -2.5000000000e-02 
               .float         -2.5000000000e-02 
C035:          .float         4.0740740741e-01  
               .float         4.0740740741e-01  
               .float         4.0740740741e-01  
               .float         4.0740740741e-01  
               .float         4.0740740741e-01  
               .float         4.0740740741e-01  
               .float         4.0740740741e-01  
               .float         4.0740740741e-01  
C036:          .float         -3.0303030303e-02 
               .float         -3.0303030303e-02 
               .float         -3.0303030303e-02 
               .float         -3.0303030303e-02 
               .float         -3.0303030303e-02 
               .float         -3.0303030303e-02 
               .float         -3.0303030303e-02 
               .float         -3.0303030303e-02 
C037:          .float         5.5000000000e-01  
               .float         5.5000000000e-01  
               .float         5.5000000000e-01  
               .float         5.5000000000e-01  
               .float         5.5000000000e-01  
               .float         5.5000000000e-01  
               .float         5.5000000000e-01  
               .float         5.5000000000e-01  
C038:          .float         -1.5625000000e-02 
               .float         -1.5625000000e-02 
               .float         -1.5625000000e-02 
               .float         -1.5625000000e-02 
               .float         -1.5625000000e-02 
               .float         -1.5625000000e-02 
               .float         -1.5625000000e-02 
               .float         -1.5625000000e-02 
C039:          .float         -1.5873015873e-02 
               .float         -1.5873015873e-02 
               .float         -1.5873015873e-02 
               .float         -1.5873015873e-02 
               .float         -1.5873015873e-02 
               .float         -1.5873015873e-02 
               .float         -1.5873015873e-02 
               .float         -1.5873015873e-02 
C040:          .float         2.7083333333e-01  
               .float         2.7083333333e-01  
               .float         2.7083333333e-01  
               .float         2.7083333333e-01  
               .float         2.7083333333e-01  
               .float         2.7083333333e-01  
               .float         2.7083333333e-01  
               .float         2.7083333333e-01  
C041:          .float         -1.6666666667e-02 
               .float         -1.6666666667e-02 
               .float         -1.6666666667e-02 
               .float         -1.6666666667e-02 
               .float         -1.6666666667e-02 
               .float         -1.6666666667e-02 
               .float         -1.6666666667e-02 
               .float         -1.6666666667e-02 
C042:          .float         2.8888888889e-01  
               .float         2.8888888889e-01  
               .float         2.8888888889e-01  
               .float         2.8888888889e-01  
               .float         2.8888888889e-01  
               .float         2.8888888889e-01  
               .float         2.8888888889e-01  
               .float         2.8888888889e-01  
C043:          .float         -1.8181818182e-02 
               .float         -1.8181818182e-02 
               .float         -1.8181818182e-02 
               .float         -1.8181818182e-02 
               .float         -1.8181818182e-02 
               .float         -1.8181818182e-02 
               .float         -1.8181818182e-02 
               .float         -1.8181818182e-02 
C044:          .float         3.2500000000e-01  
               .float         3.2500000000e-01  
               .float         3.2500000000e-01  
               .float         3.2500000000e-01  
               .float         3.2500000000e-01  
               .float         3.2500000000e-01  
               .float         3.2500000000e-01  
               .float         3.2500000000e-01  
C045:          .float         3.9393939394e-01  
               .float         3.9393939394e-01  
               .float         3.9393939394e-01  
               .float         3.9393939394e-01  
               .float         3.9393939394e-01  
               .float         3.9393939394e-01  
               .float         3.9393939394e-01  
               .float         3.9393939394e-01  
C046:          .float         -2.5641025641e-02 
               .float         -2.5641025641e-02 
               .float         -2.5641025641e-02 
               .float         -2.5641025641e-02 
               .float         -2.5641025641e-02 
               .float         -2.5641025641e-02 
               .float         -2.5641025641e-02 
               .float         -2.5641025641e-02 
C047:          .float         5.4166666667e-01  
               .float         5.4166666667e-01  
               .float         5.4166666667e-01  
               .float         5.4166666667e-01  
               .float         5.4166666667e-01  
               .float         5.4166666667e-01  
               .float         5.4166666667e-01  
               .float         5.4166666667e-01  
C048:          .float         -3.5714285714e-02 
               .float         -3.5714285714e-02 
               .float         -3.5714285714e-02 
               .float         -3.5714285714e-02 
               .float         -3.5714285714e-02 
               .float         -3.5714285714e-02 
               .float         -3.5714285714e-02 
               .float         -3.5714285714e-02 
C049:          .float         2.3809523810e-01  
               .float         2.3809523810e-01  
               .float         2.3809523810e-01  
               .float         2.3809523810e-01  
               .float         2.3809523810e-01  
               .float         2.3809523810e-01  
               .float         2.3809523810e-01  
               .float         2.3809523810e-01  
C050:          .float         2.7272727273e-01  
               .float         2.7272727273e-01  
               .float         2.7272727273e-01  
               .float         2.7272727273e-01  
               .float         2.7272727273e-01  
               .float         2.7272727273e-01  
               .float         2.7272727273e-01  
               .float         2.7272727273e-01  
C051:          .float         3.1250000000e-01  
               .float         3.1250000000e-01  
               .float         3.1250000000e-01  
               .float         3.1250000000e-01  
               .float         3.1250000000e-01  
               .float         3.1250000000e-01  
               .float         3.1250000000e-01  
               .float         3.1250000000e-01  
C052:          .float         3.8461538462e-01  
               .float         3.8461538462e-01  
               .float         3.8461538462e-01  
               .float         3.8461538462e-01  
               .float         3.8461538462e-01  
               .float         3.8461538462e-01  
               .float         3.8461538462e-01  
               .float         3.8461538462e-01  
C053:          .float         5.3571428571e-01  
               .float         5.3571428571e-01  
               .float         5.3571428571e-01  
               .float         5.3571428571e-01  
               .float         5.3571428571e-01  
               .float         5.3571428571e-01  
               .float         5.3571428571e-01  
               .float         5.3571428571e-01  
