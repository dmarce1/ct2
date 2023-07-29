
              .global        scramble

              .text

scramble:     dec            %rsi
              xor            %r8, %r8
              xor            %rdx, %rdx
main_loop:    cmp            %r8, %rdx
              jg             skip_swap
              mov            %r8, %r10
              mov            %rdx, %r11
              shl            $3, %r10
              shl            $3, %r11
              vmovaps        (%rdi, %r10, 4), %ymm0
              vmovaps        (%rdi, %r11, 4), %ymm1
              vmovaps        %ymm0, (%rdi, %r11, 4)
              vmovaps        %ymm1, (%rdi, %r10, 4)
skip_swap:    mov            %rdx, %rax
              not            %rax
              and            %rsi, %rax
              bsr            %rax, %rcx
              mov            $1, %rax
              shl            %rcx, %rax
              dec            %rax
              mov            %rdx, %rcx
              and            %rax, %rcx
              inc            %rax
              or             %rax, %rcx
              mov            %rcx, %rdx
              inc            %r8
              cmp            %r8, %rsi
              jg             main_loop
              ret

