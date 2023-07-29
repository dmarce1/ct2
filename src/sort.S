#define       ptr            %r15
#define       N              %r14
#define       lo             %r13
#define       hi             %r12
#define       xptr           %r11
#define       xmid           %r10d

              .global        particle_sort


              .text

particle_sort:push           %r15
              push           %r14
              push           %r13
              push           %r12
              push           %rbx
              mov            %rdi, ptr
              mov            %rsi, N
              mov            %rdx, lo
              mov            %rcx, hi
              mov            %r9, %rax
              mov            %r8d, xmid
              imul           N, %rax
              lea            (ptr, %rax, 4), %rax
              lea            (%rax, N, 4), %rbx
              lea            (%rbx, N, 4), %rcx
              lea            (%rcx, N, 4), %rdx
              lea            (%rdx, N, 4), %rsi
              lea            (%rsi, N, 4), %rdi
              xor            %r8, %r8
begin_loop0:  cmp            lo, hi
              jle            end_loop0
              mov            (xptr, lo, 4), %r8d
              cmp            %r8d, xmid
              jl             end_loop1
begin_loop1:  cmp            hi, lo
              je             end_loop1
              dec            hi
              mov            (xptr, hi, 4), %r8d
              cmp            %r8d, xmid
              jge            noswap
              mov            (ptr, lo, 4), %r8d
              mov            (%rax, lo, 4), %r9d
              xchg           %r8d, (ptr, hi, 4)
              xchg           %r9d, (%rax, hi, 4)
              mov            %r8d, (ptr, lo, 4)
              mov            %r9d, (%rax, lo, 4)
              mov            (%rbx, lo, 4), %r8d
              mov            (%rcx, lo, 4), %r9d
              xchg           %r8d, (%rbx, hi, 4)
              xchg           %r9d, (%rcx, hi, 4)
              mov            %r8d, (%rbx, lo, 4)
              mov            %r9d, (%rcx, lo, 4)
              mov            (%rdx, lo, 4), %r8d
              mov            (%rsi, lo, 4), %r9d
              xchg           %r8d, (%rdx, hi, 4)
              xchg           %r9d, (%rsi, hi, 4)
              mov            %r8d, (%rdx, lo, 4)
              mov            %r9d, (%rsi, lo, 4)
              mov            (%rdi, lo), %r8b
              xchg           %r8b, (%rdi, hi)
              mov            %r8b, (%rdi, lo)
noswap:       jmp            begin_loop1
end_loop1:    inc            lo
              jmp            begin_loop0
end_loop0:    mov            hi, %rax
              pop            %rbx
              pop            %r12
              pop            %r13
              pop            %r14
              pop            %r15
              ret

