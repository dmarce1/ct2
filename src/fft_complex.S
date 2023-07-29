#define        X              %r8
#define        Y              %r9
#define        N              %r10
#define        Wr             %r11
#define        Wi             %r12
#define        k2             %r13
#define        N2             %r14
#define        NTW            %r15
#define        vindex         %rsi
#define        er0            %ymm0
#define        ei0            %ymm1
#define        er1            %ymm2
#define        ei1            %ymm3
#define        er2            %ymm4
#define        ei2            %ymm5
#define        er3            %ymm6
#define        ei3            %ymm7
#define        cos1           %ymm8
#define        sin1           %ymm9
#define        cos2           %ymm10
#define        sin2           %ymm11
#define        tr             %ymm12
#define        ti             %ymm13
#define        two            %ymm14
#define        STACK_SIZE     $8
#define        vlen           -8(%rbp)

               .global        fft_complex

               .text

fft_complex:   push           %rbp
               mov            %rsp, %rbp
               sub            STACK_SIZE, %rsp
               push           %rbx
               push           %r12
               push           %r13
               push           %r14
               push           %r15
               mov            %rdi, X
               mov            %rsi, Y
               mov            %rdx, N
               mov            %rcx, vlen
               mov            N, %rdi
               push           %r8
               push           %r9
               push           %r10
               call           get_twiddles
               pop            %r10
               pop            %r9
               pop            %r8
               mov            %rax, Wr
               mov            %rdx, Wi
               mov            $1, NTW
               vmovaps        TWO, two
               call           next_level
               pop            %r15
               pop            %r14
               pop            %r13
               pop            %r12
               pop            %rbx
               mov            %rbp, %rsp
               pop            %rbp
               ret

next_level:    bsr            N, %rax
               test           $1, %rax
               jnz            do2
               call           dispatch4
               call           radix4
               jmp            done
do2:           call           dispatch2
               call           radix2
done:          ret

radix2:        mov            N, N2
               shr            N2
               mov            Y, %rdi
               sub            X, %rdi
               xor            k2, k2
oloop2:        mov            k2, %rax
               imul           NTW, %rax
               vbroadcastss   (Wr, %rax, 4), cos1
               vbroadcastss   (Wi, %rax, 4), sin1
               xor            vindex, vindex
iloop2:        mov            k2, %rax
               imul           vlen, %rax
               add            vindex, %rax
               lea            (X, %rax, 4), %rax
               lea            (%rax, N2, 4), %rbx
               vmovaps        (%rax), er0
               vmovaps        (%rbx), er1
               vmovaps        (%rax, %rdi), ei0
               vmovaps        (%rbx, %rdi), ei1
               vmovaps        er0, tr
               vmovaps        ei0, ti
               vfmadd231ps    sin1, ei1, tr
               vfnmadd231ps   sin1, er1, ti
               vfnmadd132ps   cos1, tr, er1
               vfnmadd132ps   cos1, ti, ei1
               vfmsub132ps    two, er1, er0
               vfmsub132ps    two, ei1, ei0
               vmovaps        er0, (%rax)
               vmovaps        er1, (%rbx)
               vmovaps        ei0, (%rax, %rdi)
               vmovaps        ei1, (%rbx, %rdi)
               add            $8, vindex
               cmp            vlen, vindex
               jne            iloop2
               inc            k2
               cmp            k2, N2
               jne            oloop2
               ret

radix4:        mov            N, N2
               shr            $2, N2
               mov            Y, %rdi
               sub            X, %rdi
               xor            k2, k2
oloop4:        mov            k2, %rax
               imul           NTW, %rax
               imul           $2, %rax, %rdx
               vbroadcastss   (Wr, %rax, 4), cos1
               vbroadcastss   (Wi, %rax, 4), sin1
               vbroadcastss   (Wr, %rdx, 4), cos2
               vbroadcastss   (Wi, %rdx, 4), sin2
               xor            vindex, vindex
iloop4:        mov            k2, %rax
               imul           vlen, %rax
               add            vindex, %rax
               lea            (X, %rax, 4), %rax
               lea            (%rax, N2, 4), %rbx
               lea            (%rbx, N2, 4), %rcx
               lea            (%rcx, N2, 4), %rdx
               vmovaps        (%rax), er0
               vmovaps        (%rbx), er1
               vmovaps        (%rcx), er2
               vmovaps        (%rdx), er3
               vmovaps        (%rax, %rdi), ei0
               vmovaps        (%rbx, %rdi), ei1
               vmovaps        (%rcx, %rdi), ei2
               vmovaps        (%rdx, %rdi), ei3
               vmovaps        er0, tr
               vmovaps        ei0, ti
               vfmadd231ps    sin2, ei2, tr
               vfnmadd231ps   sin2, er2, ti
               vfnmadd132ps   cos2, tr, er2
               vfnmadd132ps   cos2, ti, ei2
               vfmsub132ps    two, er2, er0
               vfmsub132ps    two, ei2, ei0
               vmovaps        er1, tr
               vmovaps        ei1, ti
               vfmadd231ps    sin2, ei3, tr
               vfnmadd231ps   sin2, er3, ti
               vfnmadd132ps   cos2, tr, er3
               vfnmadd132ps   cos2, ti, ei3
               vfmsub132ps    two, er3, er1
               vfmsub132ps    two, ei3, ei1
               vmovaps        er0, tr
               vmovaps        ei0, ti
               vfmadd231ps    sin1, ei1, tr
               vfnmadd231ps   sin1, er1, ti
               vfnmadd132ps   cos1, tr, er1
               vfmsub132ps    cos1, ti, ei1
               vfmsub132ps    two, er1, er0
               vfmadd132ps    two, ei1, ei0
               vmovaps        er2, tr
               vmovaps        ei2, ti
               vfmadd231ps    cos1, ei3, tr
               vfnmadd231ps   cos1, er3, ti
               vfmadd132ps    sin1, tr, er3
               vfmadd132ps    sin1, ti, ei3
               vfmsub132ps    two, er3, er2
               vfnmadd132ps   two, ei3, ei2
               vmovaps        er0, (%rax)
               vmovaps        er3, (%rbx)
               vmovaps        er1, (%rcx)
               vmovaps        er2, (%rdx)
               vmovaps        ei0, (%rax, %rdi)
               vmovaps        ei3, (%rbx, %rdi)
               vmovaps        ei1, (%rcx, %rdi)
               vmovaps        ei2, (%rdx, %rdi)
               add            $8, vindex
               cmp            vlen, vindex
               jne            iloop4
               inc            k2
               cmp            k2, N2
               jne            oloop4
               ret

dispatch2:     cmp            $2, N
               jle            no2
               mov            N, %rax
               shr            %rax
               imul           vlen, %rax
               push           X
               push           Y
               lea            (X, %rax, 4), X
               lea            (Y, %rax, 4), Y
               shr            N
               shl            NTW
               call           next_level
               pop            Y
               pop            X
               call           next_level
               shl            N
               shr            NTW
no2:           ret

dispatch4:     cmp            $4, N
               jle            no4
               mov            N, %rax
               shr            $2, %rax
               imul           vlen, %rax
               push           X
               push           Y
               lea            (X, %rax, 4), X
               lea            (Y, %rax, 4), Y
               push           X
               push           Y
               lea            (X, %rax, 4), X
               lea            (Y, %rax, 4), Y
               push           X
               push           Y
               lea            (X, %rax, 4), X
               lea            (Y, %rax, 4), Y
               shr            $2, N
               shl            $2, NTW
               call           next_level
               pop            Y
               pop            X
               call           next_level
               pop            Y
               pop            X
               call           next_level
               pop            Y
               pop            X
               call           next_level
               shl            $2, N
               shr            $2, NTW
no4:           ret

               .align         32
TWO:           .float         2.0
               .float         2.0
               .float         2.0
               .float         2.0
               .float         2.0
               .float         2.0
               .float         2.0
               .float         2.0













