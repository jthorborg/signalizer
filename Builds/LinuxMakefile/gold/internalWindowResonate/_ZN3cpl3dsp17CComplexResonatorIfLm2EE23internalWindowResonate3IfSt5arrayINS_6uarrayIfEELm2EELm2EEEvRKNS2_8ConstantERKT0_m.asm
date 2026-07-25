; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000964940 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long):
  964940:	endbr64
  964944:	push   r15
  964946:	push   r14
  964948:	push   r13
  96494a:	mov    r13,rsi
  96494d:	push   r12
  96494f:	mov    r12,rdi
  964952:	push   rbp
  964953:	push   rbx
  964954:	mov    rbx,rcx
  964957:	sub    rsp,0x58
  96495b:	mov    QWORD PTR [rsp],rdx
  96495f:	movzx  ebp,BYTE PTR [rip+0x2cd76d]        # c320d3 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  964966:	test   bpl,bpl
  964969:	je     964de9 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4a9>
  96496f:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  96497b:	cmp    QWORD PTR [rax-0x1d8],0x0
  964983:	mov    r14,rax
  964986:	lea    rcx,[rax-0x360]
  96498d:	je     9649b8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x78>
  96498f:	movzx  r15d,BYTE PTR [rax-0x1e0]
  964997:	cmp    r15b,0xf
  96499b:	jbe    964e4f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x50f>
  9649a1:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  9649ad:	add    r15d,0x1
  9649b1:	mov    BYTE PTR [rax-0x1e0],r15b
  9649b8:	mov    rdx,QWORD PTR [r13+0x38]
  9649bc:	mov    r8,QWORD PTR [r13+0x48]
  9649c0:	mov    r9,QWORD PTR [r12+0x8]
  9649c5:	mov    rsi,QWORD PTR [r12]
  9649c9:	mov    rdi,rdx
  9649cc:	imul   rdi,r8
  9649d0:	mov    rax,r9
  9649d3:	sub    rax,rsi
  9649d6:	sar    rax,0x2
  9649da:	lea    rcx,[rdi*4+0x0]
  9649e2:	cmp    rax,rcx
  9649e5:	jb     964e34 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4f4>
  9649eb:	cmp    rcx,rax
  9649ee:	jb     964dcf <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x48f>
  9649f4:	mov    rsi,QWORD PTR [r13+0x40]
  9649f8:	lea    rax,[r8+r8*1]
  9649fc:	imul   rdx,rax
  964a00:	mov    QWORD PTR [rsp+0x28],rsi
  964a05:	test   rsi,rsi
  964a08:	je     964ca0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x360>
  964a0e:	mov    rsi,QWORD PTR [rsp]
  964a12:	lea    rcx,[r8*4+0x0]
  964a1a:	shl    r8,0x3
  964a1e:	mov    r9,QWORD PTR [r12]
  964a22:	lea    r11,[rcx+rdx*1]
  964a26:	mov    r14,QWORD PTR [r13+0x0]
  964a2a:	mov    rdi,QWORD PTR [rsi]
  964a2d:	mov    rsi,QWORD PTR [rsi+0x10]
  964a31:	lea    r12,[r9+rdx*4]
  964a35:	mov    rdx,r11
  964a38:	lea    r13,[r14+r8*1]
  964a3c:	sub    rdx,rax
  964a3f:	lea    r10,[r9+r8*1]
  964a43:	lea    r15,[r13+r8*1+0x0]
  964a48:	add    r8,r10
  964a4b:	lea    rbp,[r9+rdx*4]
  964a4f:	xor    edx,edx
  964a51:	lea    r11,[r9+r11*4]
  964a55:	nop    DWORD PTR [rax]
  964a58:	vmovss xmm1,DWORD PTR [r9+rcx*1]
  964a5e:	vmovss xmm0,DWORD PTR [r8+rcx*1]
  964a64:	vmovss xmm2,DWORD PTR [r10+rcx*1]
  964a6a:	vmovss DWORD PTR [rsp+0x24],xmm1
  964a70:	vmovss xmm12,DWORD PTR [r11+rcx*1]
  964a76:	vmovss DWORD PTR [rsp+0x18],xmm0
  964a7c:	vmovss xmm1,DWORD PTR [r12+rcx*1]
  964a82:	vmovss DWORD PTR [rsp],xmm2
  964a87:	vmovss xmm0,DWORD PTR [rbp+rcx*1+0x0]
  964a8d:	vmovss DWORD PTR [rsp+0x10],xmm12
  964a93:	vmovss DWORD PTR [rsp+0x20],xmm1
  964a99:	vmovss xmm6,DWORD PTR [r14+rdx*4]
  964a9f:	vmovss DWORD PTR [rsp+0x1c],xmm0
  964aa5:	vmovss xmm7,DWORD PTR [r14+rcx*1]
  964aab:	vmovss xmm9,DWORD PTR [r13+rdx*4+0x0]
  964ab2:	vmovss xmm10,DWORD PTR [r13+rcx*1+0x0]
  964ab9:	vmovss xmm11,DWORD PTR [r15+rdx*4]
  964abf:	vmovss xmm8,DWORD PTR [r15+rcx*1]
  964ac5:	vmovss xmm5,DWORD PTR [r9+rdx*4]
  964acb:	vmovss xmm4,DWORD PTR [r10+rdx*4]
  964ad1:	vmovss xmm3,DWORD PTR [r8+rdx*4]
  964ad7:	vmovss xmm2,DWORD PTR [r12+rdx*4]
  964add:	vmovss xmm1,DWORD PTR [rbp+rdx*4+0x0]
  964ae3:	vmovss xmm0,DWORD PTR [r11+rdx*4]
  964ae9:	test   rbx,rbx
  964aec:	je     964c1f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x2df>
  964af2:	xor    eax,eax
  964af4:	nop    DWORD PTR [rax+0x0]
  964af8:	vmovss xmm15,DWORD PTR [rsp+0x24]
  964afe:	vmulss xmm14,xmm5,xmm7
  964b02:	vmulss xmm5,xmm5,xmm6
  964b06:	vmovss xmm12,DWORD PTR [rdi+rax*4]
  964b0b:	vmulss xmm13,xmm7,xmm15
  964b10:	vmulss xmm15,xmm6,xmm15
  964b15:	vaddss xmm5,xmm5,xmm12
  964b1a:	vaddss xmm15,xmm14,xmm15
  964b1f:	vmovss xmm14,DWORD PTR [rsp]
  964b24:	vsubss xmm5,xmm5,xmm13
  964b29:	vmulss xmm13,xmm10,xmm14
  964b2e:	vmovss DWORD PTR [rsp+0x24],xmm15
  964b34:	vmulss xmm14,xmm9,xmm14
  964b39:	vmulss xmm15,xmm10,xmm4
  964b3d:	vmulss xmm4,xmm9,xmm4
  964b41:	vaddss xmm15,xmm15,xmm14
  964b46:	vmovss xmm14,DWORD PTR [rsp+0x18]
  964b4c:	vaddss xmm4,xmm4,xmm12
  964b51:	vmovss DWORD PTR [rsp],xmm15
  964b56:	vmulss xmm15,xmm8,xmm3
  964b5a:	vsubss xmm4,xmm4,xmm13
  964b5f:	vmulss xmm3,xmm11,xmm3
  964b63:	vmulss xmm13,xmm8,xmm14
  964b68:	vmulss xmm14,xmm11,xmm14
  964b6d:	vaddss xmm3,xmm3,xmm12
  964b72:	vmovss xmm12,DWORD PTR [rsi+rax*4]
  964b77:	add    rax,0x1
  964b7b:	vaddss xmm15,xmm15,xmm14
  964b80:	vmovss xmm14,DWORD PTR [rsp+0x20]
  964b86:	vsubss xmm3,xmm3,xmm13
  964b8b:	vmulss xmm13,xmm7,xmm14
  964b90:	vmovss DWORD PTR [rsp+0x18],xmm15
  964b96:	vmulss xmm14,xmm6,xmm14
  964b9b:	vmulss xmm15,xmm7,xmm2
  964b9f:	vmulss xmm2,xmm6,xmm2
  964ba3:	vaddss xmm15,xmm15,xmm14
  964ba8:	vmovss xmm14,DWORD PTR [rsp+0x1c]
  964bae:	vaddss xmm2,xmm2,xmm12
  964bb3:	vmovss DWORD PTR [rsp+0x20],xmm15
  964bb9:	vmulss xmm15,xmm1,xmm10
  964bbe:	vsubss xmm2,xmm2,xmm13
  964bc3:	vmulss xmm1,xmm1,xmm9
  964bc8:	vmulss xmm13,xmm10,xmm14
  964bcd:	vmulss xmm14,xmm9,xmm14
  964bd2:	vaddss xmm1,xmm1,xmm12
  964bd7:	vaddss xmm15,xmm15,xmm14
  964bdc:	vmovss xmm14,DWORD PTR [rsp+0x10]
  964be2:	vsubss xmm1,xmm1,xmm13
  964be7:	vmulss xmm13,xmm8,xmm14
  964bec:	vmovss DWORD PTR [rsp+0x1c],xmm15
  964bf2:	vmulss xmm14,xmm11,xmm14
  964bf7:	vmulss xmm15,xmm0,xmm8
  964bfc:	vmulss xmm0,xmm0,xmm11
  964c01:	vaddss xmm15,xmm15,xmm14
  964c06:	vaddss xmm0,xmm0,xmm12
  964c0b:	vmovss DWORD PTR [rsp+0x10],xmm15
  964c11:	vsubss xmm0,xmm0,xmm13
  964c16:	cmp    rbx,rax
  964c19:	jne    964af8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x1b8>
  964c1f:	vmovss xmm6,DWORD PTR [rsp+0x24]
  964c25:	vmovss DWORD PTR [r9+rdx*4],xmm5
  964c2b:	vmovss xmm7,DWORD PTR [rsp]
  964c30:	vmovss DWORD PTR [r9+rcx*1],xmm6
  964c36:	vmovss xmm6,DWORD PTR [rsp+0x18]
  964c3c:	vmovss DWORD PTR [r10+rdx*4],xmm4
  964c42:	mov    rax,QWORD PTR [rsp+0x28]
  964c47:	vmovss DWORD PTR [r10+rcx*1],xmm7
  964c4d:	vmovss xmm7,DWORD PTR [rsp+0x20]
  964c53:	vmovss DWORD PTR [r8+rdx*4],xmm3
  964c59:	vmovss DWORD PTR [r8+rcx*1],xmm6
  964c5f:	vmovss xmm6,DWORD PTR [rsp+0x1c]
  964c65:	vmovss DWORD PTR [r12+rdx*4],xmm2
  964c6b:	vmovss DWORD PTR [r12+rcx*1],xmm7
  964c71:	vmovss xmm7,DWORD PTR [rsp+0x10]
  964c77:	vmovss DWORD PTR [rbp+rdx*4+0x0],xmm1
  964c7d:	vmovss DWORD PTR [rbp+rcx*1+0x0],xmm6
  964c83:	vmovss DWORD PTR [r11+rdx*4],xmm0
  964c89:	add    rdx,0x1
  964c8d:	vmovss DWORD PTR [r11+rcx*1],xmm7
  964c93:	add    rcx,0x4
  964c97:	cmp    rdx,rax
  964c9a:	jne    964a58 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x118>
  964ca0:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  964cac:	lea    rbp,[rax-0x360]
  964cb3:	mov    r12,rax
  964cb6:	mov    r13,QWORD PTR [rbp+0x188]
  964cbd:	test   r13,r13
  964cc0:	je     964ce0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3a0>
  964cc2:	movzx  eax,BYTE PTR [rbp+0x180]
  964cc9:	lea    ebx,[rax-0x1]
  964ccc:	mov    BYTE PTR [rbp+0x180],bl
  964cd2:	cmp    bl,0xf
  964cd5:	jbe    964cef <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3af>
  964cd7:	add    WORD PTR [r13+0xfe2],0x1
  964ce0:	add    rsp,0x58
  964ce4:	pop    rbx
  964ce5:	pop    rbp
  964ce6:	pop    r12
  964ce8:	pop    r13
  964cea:	pop    r14
  964cec:	pop    r15
  964cee:	ret
  964cef:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  964cf4:	movzx  r15d,bl
  964cf8:	lea    rdx,[r15+r15*2]
  964cfc:	mov    rcx,rax
  964cff:	vmovq  xmm6,rax
  964d04:	shl    rdx,0x3
  964d08:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  964d0d:	vpinsrq xmm0,xmm6,rcx,0x1
  964d13:	vpsubq xmm0,xmm0,XMMWORD PTR [r12+rdx*1-0x360]
  964d1d:	movzx  r12d,BYTE PTR [rbp+0x181]
  964d25:	cmp    r12b,bl
  964d28:	jae    964d3e <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3fe>
  964d2a:	lea    eax,[r15-0x1]
  964d2e:	vmovq  rdx,xmm0
  964d33:	cdqe
  964d35:	lea    rax,[rax+rax*2]
  964d39:	add    QWORD PTR [rbp+rax*8+0x8],rdx
  964d3e:	movzx  ebp,WORD PTR [r13+0xfe0]
  964d46:	cmp    bp,0x7f
  964d4a:	je     964e9b <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x55b>
  964d50:	vmovdqa XMMWORD PTR [rsp],xmm0
  964d55:	sub    ebx,r12d
  964d58:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  964d64:	lea    rdx,[r15+r15*2]
  964d68:	vmovdqa xmm0,XMMWORD PTR [rsp]
  964d6d:	mov    DWORD PTR [rsp+0x48],0x0
  964d75:	vmovdqu XMMWORD PTR [rsp+0x38],xmm0
  964d7b:	lea    rdx,[rax+rdx*8-0x360]
  964d83:	mov    rax,QWORD PTR [rdx]
  964d86:	mov    QWORD PTR [rsp+0x30],rax
  964d8b:	movzx  eax,BYTE PTR [rdx+0x10]
  964d8f:	lea    edx,[rbp+0x1]
  964d92:	vmovdqa xmm7,XMMWORD PTR [rsp+0x30]
  964d98:	mov    WORD PTR [r13+0xfe0],dx
  964da0:	mov    ah,bl
  964da2:	mov    WORD PTR [rsp+0x4c],ax
  964da7:	movzx  eax,bp
  964daa:	vmovdqu xmm6,XMMWORD PTR [rsp+0x3e]
  964db0:	shl    rax,0x5
  964db4:	add    rax,r13
  964db7:	vmovdqu XMMWORD PTR [rax],xmm7
  964dbb:	vmovdqu XMMWORD PTR [rax+0xe],xmm6
  964dc0:	add    rsp,0x58
  964dc4:	pop    rbx
  964dc5:	pop    rbp
  964dc6:	pop    r12
  964dc8:	pop    r13
  964dca:	pop    r14
  964dcc:	pop    r15
  964dce:	ret
  964dcf:	shl    rdi,0x4
  964dd3:	add    rsi,rdi
  964dd6:	cmp    r9,rsi
  964dd9:	je     9649f4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb4>
  964ddf:	mov    QWORD PTR [r12+0x8],rsi
  964de4:	jmp    9649f4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb4>
  964de9:	mov    eax,0x1
  964dee:	lock xadd DWORD PTR [rip+0x2c1d6a],eax        # c26b60 <cpl::Profiling::registerRegion(char const*)::counter>
  964df6:	mov    ecx,0x1
  964dfb:	mov    ebp,0x1
  964e00:	add    eax,0x2
  964e03:	cmp    eax,0xfe
  964e08:	ja     964e22 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4e2>
  964e0a:	lea    rdx,[rip+0x2c020f]        # c25020 <cpl::Profiling::regions>
  964e11:	mov    ecx,eax
  964e13:	mov    ebp,eax
  964e15:	lea    rdi,[rip+0x9b496]        # a002b2 <_IO_stdin_used+0x182b2>
  964e1c:	mov    QWORD PTR [rdx+rcx*8],rdi
  964e20:	mov    ecx,eax
  964e22:	xor    eax,eax
  964e24:	lock cmpxchg BYTE PTR [rip+0x2cd2a7],cl        # c320d3 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  964e2c:	cmovne ebp,eax
  964e2f:	jmp    96496f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x2f>
  964e34:	mov    rsi,rcx
  964e37:	mov    rdi,r12
  964e3a:	sub    rsi,rax
  964e3d:	call   939b50 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  964e42:	mov    rdx,QWORD PTR [r13+0x38]
  964e46:	mov    r8,QWORD PTR [r13+0x48]
  964e4a:	jmp    9649f4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb4>
  964e4f:	mov    QWORD PTR [rsp+0x10],rcx
  964e54:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  964e59:	mov    rcx,QWORD PTR [rsp+0x10]
  964e5e:	mov    BYTE PTR [rsp+0x40],bpl
  964e63:	mov    QWORD PTR [rsp+0x30],rax
  964e68:	lea    rax,[r15+r15*2]
  964e6c:	mov    QWORD PTR [rsp+0x38],0x0
  964e75:	vmovdqa xmm7,XMMWORD PTR [rsp+0x30]
  964e7b:	shl    rax,0x3
  964e7f:	vmovdqu XMMWORD PTR [r14+rax*1-0x360],xmm7
  964e89:	mov    BYTE PTR [rax+rcx*1+0x10],bpl
  964e8e:	movzx  r15d,BYTE PTR [rcx+0x180]
  964e96:	jmp    9649a1 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x61>
  964e9b:	movzx  eax,WORD PTR [r13+0xfe2]
  964ea3:	cmp    ax,0xffff
  964ea7:	je     964ce0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3a0>
  964ead:	add    eax,0x1
  964eb0:	mov    WORD PTR [r13+0xfe2],ax
  964eb8:	jmp    964ce0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3a0>
  964ebd:	endbr64
  964ec1:	mov    rbx,rax
  964ec4:	vzeroupper
  964ec7:	call   93c730 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  964ecc:	mov    rdi,rbx
  964ecf:	call   e1b10 <_Unwind_Resume@plt>

Disassembly of section .fini:
