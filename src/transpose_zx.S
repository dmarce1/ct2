#define       X              %r8
#define       M              %r9
#define       Ix             %r10
#define       Iy             %r11
#define       Iz             %r12
#define       N              %r13
#define       M0             $4

              .global        transpose_zx


              .text

transpose_zx: push           %r13
              push           %r12
              push           %rbx
              mov            %rdi, X
              mov            %rsi, N
              xor            Ix, Ix
              xor            Iy, Iy
              xor            Iz, Iz
              mov            N, M
              mov            N, %rax
              imul           %rax, %rax
              imul           $2, %rax, %rbx
              imul           $3, %rax, %rcx
              call           next_level
              pop            %rbx
              pop            %r12
              pop            %r13
              ret
next_level:   cmp            Ix, Iz
              jl             done
              cmp            M0, M
              jg             next_levels
              cmp            Ix, Iz
              je             on_diagonal
              mov            Ix, %rsi
              imul           N, %rsi
              add            Iy, %rsi
              imul           N, %rsi
              add            Iz, %rsi
              mov            Iz, %rdi
              imul           N, %rdi
              add            Iy, %rdi
              imul           N, %rdi
              add            Ix, %rdi
              lea            (X, %rsi, 4), %rsi
              lea            (X, %rdi, 4), %rdi
              vmovapd        (%rsi), %ymm0
              vmovapd        (%rsi, %rax, 4), %ymm1
              vmovapd        (%rsi, %rbx, 4), %ymm2
              vmovapd        (%rsi, %rcx, 4), %ymm3
              vmovapd        (%rdi), %ymm8
              vmovapd        (%rdi, %rax, 4), %ymm9
              vmovapd        (%rdi, %rbx, 4), %ymm10
              vmovapd        (%rdi, %rcx, 4), %ymm11
              vunpcklpd      %ymm1, %ymm0, %ymm4
              vunpckhpd      %ymm1, %ymm0, %ymm5
              vunpcklpd      %ymm3, %ymm2, %ymm6
              vunpckhpd      %ymm3, %ymm2, %ymm7
              vperm2f128     $0x20, %ymm6, %ymm4, %ymm0
              vperm2f128     $0x20, %ymm7, %ymm5, %ymm1
              vperm2f128     $0x31, %ymm6, %ymm4, %ymm2
              vperm2f128     $0x31, %ymm7, %ymm5, %ymm3
              vunpcklpd      %ymm9, %ymm8, %ymm4
              vunpckhpd      %ymm9, %ymm8, %ymm5
              vunpcklpd      %ymm11, %ymm10, %ymm6
              vunpckhpd      %ymm11, %ymm10, %ymm7
              vperm2f128     $0x20, %ymm6, %ymm4, %ymm8
              vperm2f128     $0x20, %ymm7, %ymm5, %ymm9
              vperm2f128     $0x31, %ymm6, %ymm4, %ymm10
              vperm2f128     $0x31, %ymm7, %ymm5, %ymm11
              vmovapd        %ymm0, (%rdi)
              vmovapd        %ymm1, (%rdi, %rax, 4)
              vmovapd        %ymm2, (%rdi, %rbx, 4)
              vmovapd        %ymm3, (%rdi, %rcx, 4)
              vmovapd        %ymm8, (%rsi)
              vmovapd        %ymm9, (%rsi, %rax, 4)
              vmovapd        %ymm10, (%rsi, %rbx, 4)
              vmovapd        %ymm11, (%rsi, %rcx, 4)
              jmp            done
on_diagonal:  mov            Ix, %rsi
              imul           N, %rsi
              add            Iy, %rsi
              imul           N, %rsi
              add            Iz, %rsi
              lea            (X, %rsi, 4), %rsi
              vmovapd        (%rsi), %ymm0
              vmovapd        (%rsi, %rax, 4), %ymm1
              vmovapd        (%rsi, %rbx, 4), %ymm2
              vmovapd        (%rsi, %rcx, 4), %ymm3
              vunpcklpd      %ymm1, %ymm0, %ymm4
              vunpckhpd      %ymm1, %ymm0, %ymm5
              vunpcklpd      %ymm3, %ymm2, %ymm6
              vunpckhpd      %ymm3, %ymm2, %ymm7
              vperm2f128     $0x20, %ymm6, %ymm4, %ymm0
              vperm2f128     $0x20, %ymm7, %ymm5, %ymm1
              vperm2f128     $0x31, %ymm6, %ymm4, %ymm2
              vperm2f128     $0x31, %ymm7, %ymm5, %ymm3
              vmovapd        %ymm0, (%rsi)
              vmovapd        %ymm1, (%rsi, %rax, 4)
              vmovapd        %ymm2, (%rsi, %rbx, 4)
              vmovapd        %ymm3, (%rsi, %rcx, 4)
              jmp            done
next_levels:  shr            M
              call           next_level
              add            M, Iy
              call           next_level
              sub            M, Iy
              add            M, Ix
              call           next_level
              add            M, Iy
              call           next_level
              sub            M, Iy
              sub            M, Ix
              add            M, Iz
              call           next_level
              add            M, Iy
              call           next_level
              sub            M, Iy
              add            M, Ix
              call           next_level
              add            M, Iy
              call           next_level
              sub            M, Iy
              sub            M, Ix
              sub            M, Iz
              shl            M
done:         ret





