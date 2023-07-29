#define       logN           %r12
#define       LOCK           %r13
#define       C              %r14
#define       S              %r15
#define       STACK_SIZE     $64
#define       N2m1           -8(%rbp)
#define       N              -16(%rbp)
#define       twopioN        -24(%rbp)
#define       cos            -32(%rbp)
#define       sin            -40(%rbp)
#define       N2             -48(%rbp)

              .global        get_twiddles


              .text

get_twiddles: push           %rbp
              mov            %rsp, %rbp
              sub            STACK_SIZE, %rsp
              push           %rbx
              push           %r15
              push           %r12
              push           %r13
              push           %r14
              mov            %rdi, N
              bsr            N, logN
              mov            N, %rax
              mov            %rax, N2
              dec            %rax
              mov            %rax, N2m1
              lea            (locks), %rax
              lea            (%rax, logN, 8), LOCK
spin_loop:    xor            %rax, %rax
              xchg           (LOCK), %rax
              cmp            $1, %rax
              je             lock_acquired
              pause
              jmp            spin_loop
lock_acquired:lea            (Cs), %rax
              lea            (%rax, logN, 8), %rax
              mov            (%rax), %rax
              cmp            $0, %rax
              jne            load_from_mem
              imul           $4, N, %rdi
              add            $32, %rdi
              call           malloc
              add            $32, %rax
              and            $0xffffffffffffffe0, %rax
              mov            %rax, C
              imul           $4, N, %rdi
              add            $32, %rdi
              call           malloc
              add            $32, %rax
              and            $0xffffffffffffffe0, %rax
              mov            %rax, S
              emms
              fldpi
              fldpi
              faddp
              fildl          N
              fdivrp
              fchs
              fstpl          twopioN
              xor            %rcx, %rcx
main_loop:    mov            %rcx, cos
              fildl          cos
              fmull          twopioN
              fcos
              fstps          (C, %rcx, 4)
              mov            %rcx, sin
              fildl          sin
              fmull          twopioN
              fsin
              fstps          (S, %rcx, 4)
              cmp            %rcx, N2m1
              je             load_complete
              inc            %rcx
              jmp            main_loop
load_complete:lea            (Cs), %rax
              lea            (%rax, logN, 8), %rax
              mov            C, (%rax)
              lea            (Ss), %rax
              lea            (%rax, logN, 8), %rax
              mov            S, (%rax)
              jmp            done
load_from_mem:lea            (Cs), %rax
              lea            (%rax, logN, 8), %rax
              mov            (%rax), C
              lea            (Ss), %rax
              lea            (%rax, logN, 8), %rax
              mov            (%rax), S
done:         mov            C, %rax
              mov            S, %rdx
              mov            $1, %rcx
              xchg           (LOCK), %rcx
              pop            %r14
              pop            %r13
              pop            %r12
              pop            %r15
              pop            %rbx
              mov            %rbp, %rsp
              pop            %rbp
              ret

              .data

              .align         8
Cs:           .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
Ss:           .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
              .quad          0
locks:        .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
              .quad          1
