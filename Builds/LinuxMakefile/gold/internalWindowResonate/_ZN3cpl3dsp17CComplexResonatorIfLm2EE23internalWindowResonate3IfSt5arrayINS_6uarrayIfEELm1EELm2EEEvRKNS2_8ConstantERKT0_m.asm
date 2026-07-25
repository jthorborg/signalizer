; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000960c20 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long):
  960c20:	endbr64
  960c24:	push   r15
  960c26:	push   r14
  960c28:	push   r13
  960c2a:	mov    r13,rsi
  960c2d:	push   r12
  960c2f:	mov    r12,rdi
  960c32:	push   rbp
  960c33:	push   rbx
  960c34:	mov    rbx,rcx
  960c37:	sub    rsp,0x58
  960c3b:	mov    QWORD PTR [rsp],rdx
  960c3f:	movzx  ebp,BYTE PTR [rip+0x2d1497]        # c320dd <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  960c46:	test   bpl,bpl
  960c49:	je     9610c9 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4a9>
  960c4f:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  960c5b:	cmp    QWORD PTR [rax-0x1d8],0x0
  960c63:	mov    r14,rax
  960c66:	lea    rcx,[rax-0x360]
  960c6d:	je     960c98 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x78>
  960c6f:	movzx  r15d,BYTE PTR [rax-0x1e0]
  960c77:	cmp    r15b,0xf
  960c7b:	jbe    96112f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x50f>
  960c81:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  960c8d:	add    r15d,0x1
  960c91:	mov    BYTE PTR [rax-0x1e0],r15b
  960c98:	mov    rdx,QWORD PTR [r13+0x38]
  960c9c:	mov    r8,QWORD PTR [r13+0x48]
  960ca0:	mov    r9,QWORD PTR [r12+0x8]
  960ca5:	mov    rsi,QWORD PTR [r12]
  960ca9:	mov    rdi,rdx
  960cac:	imul   rdi,r8
  960cb0:	mov    rax,r9
  960cb3:	sub    rax,rsi
  960cb6:	sar    rax,0x2
  960cba:	lea    rcx,[rdi*4+0x0]
  960cc2:	cmp    rax,rcx
  960cc5:	jb     961114 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4f4>
  960ccb:	cmp    rcx,rax
  960cce:	jb     9610af <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x48f>
  960cd4:	mov    rsi,QWORD PTR [r13+0x40]
  960cd8:	lea    rax,[r8+r8*1]
  960cdc:	imul   rdx,rax
  960ce0:	mov    QWORD PTR [rsp+0x28],rsi
  960ce5:	test   rsi,rsi
  960ce8:	je     960f80 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x360>
  960cee:	mov    rsi,QWORD PTR [rsp]
  960cf2:	lea    rcx,[r8*4+0x0]
  960cfa:	shl    r8,0x3
  960cfe:	mov    r9,QWORD PTR [r12]
  960d02:	lea    r11,[rcx+rdx*1]
  960d06:	mov    r14,QWORD PTR [r13+0x0]
  960d0a:	mov    rdi,QWORD PTR [rsi]
  960d0d:	mov    rsi,QWORD PTR [rsi+0x10]
  960d11:	lea    r12,[r9+rdx*4]
  960d15:	mov    rdx,r11
  960d18:	lea    r13,[r14+r8*1]
  960d1c:	sub    rdx,rax
  960d1f:	lea    r10,[r9+r8*1]
  960d23:	lea    r15,[r13+r8*1+0x0]
  960d28:	add    r8,r10
  960d2b:	lea    rbp,[r9+rdx*4]
  960d2f:	xor    edx,edx
  960d31:	lea    r11,[r9+r11*4]
  960d35:	nop    DWORD PTR [rax]
  960d38:	vmovss xmm1,DWORD PTR [r9+rcx*1]
  960d3e:	vmovss xmm0,DWORD PTR [r8+rcx*1]
  960d44:	vmovss xmm2,DWORD PTR [r10+rcx*1]
  960d4a:	vmovss DWORD PTR [rsp+0x24],xmm1
  960d50:	vmovss xmm12,DWORD PTR [r11+rcx*1]
  960d56:	vmovss DWORD PTR [rsp+0x18],xmm0
  960d5c:	vmovss xmm1,DWORD PTR [r12+rcx*1]
  960d62:	vmovss DWORD PTR [rsp],xmm2
  960d67:	vmovss xmm0,DWORD PTR [rbp+rcx*1+0x0]
  960d6d:	vmovss DWORD PTR [rsp+0x10],xmm12
  960d73:	vmovss DWORD PTR [rsp+0x20],xmm1
  960d79:	vmovss xmm6,DWORD PTR [r14+rdx*4]
  960d7f:	vmovss DWORD PTR [rsp+0x1c],xmm0
  960d85:	vmovss xmm7,DWORD PTR [r14+rcx*1]
  960d8b:	vmovss xmm9,DWORD PTR [r13+rdx*4+0x0]
  960d92:	vmovss xmm10,DWORD PTR [r13+rcx*1+0x0]
  960d99:	vmovss xmm11,DWORD PTR [r15+rdx*4]
  960d9f:	vmovss xmm8,DWORD PTR [r15+rcx*1]
  960da5:	vmovss xmm5,DWORD PTR [r9+rdx*4]
  960dab:	vmovss xmm4,DWORD PTR [r10+rdx*4]
  960db1:	vmovss xmm3,DWORD PTR [r8+rdx*4]
  960db7:	vmovss xmm2,DWORD PTR [r12+rdx*4]
  960dbd:	vmovss xmm1,DWORD PTR [rbp+rdx*4+0x0]
  960dc3:	vmovss xmm0,DWORD PTR [r11+rdx*4]
  960dc9:	test   rbx,rbx
  960dcc:	je     960eff <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x2df>
  960dd2:	xor    eax,eax
  960dd4:	nop    DWORD PTR [rax+0x0]
  960dd8:	vmovss xmm15,DWORD PTR [rsp+0x24]
  960dde:	vmulss xmm14,xmm5,xmm7
  960de2:	vmulss xmm5,xmm5,xmm6
  960de6:	vmovss xmm12,DWORD PTR [rdi+rax*4]
  960deb:	vmulss xmm13,xmm7,xmm15
  960df0:	vmulss xmm15,xmm6,xmm15
  960df5:	vaddss xmm5,xmm5,xmm12
  960dfa:	vaddss xmm15,xmm14,xmm15
  960dff:	vmovss xmm14,DWORD PTR [rsp]
  960e04:	vsubss xmm5,xmm5,xmm13
  960e09:	vmulss xmm13,xmm10,xmm14
  960e0e:	vmovss DWORD PTR [rsp+0x24],xmm15
  960e14:	vmulss xmm14,xmm9,xmm14
  960e19:	vmulss xmm15,xmm10,xmm4
  960e1d:	vmulss xmm4,xmm9,xmm4
  960e21:	vaddss xmm15,xmm15,xmm14
  960e26:	vmovss xmm14,DWORD PTR [rsp+0x18]
  960e2c:	vaddss xmm4,xmm4,xmm12
  960e31:	vmovss DWORD PTR [rsp],xmm15
  960e36:	vmulss xmm15,xmm8,xmm3
  960e3a:	vsubss xmm4,xmm4,xmm13
  960e3f:	vmulss xmm3,xmm11,xmm3
  960e43:	vmulss xmm13,xmm8,xmm14
  960e48:	vmulss xmm14,xmm11,xmm14
  960e4d:	vaddss xmm3,xmm3,xmm12
  960e52:	vmovss xmm12,DWORD PTR [rsi+rax*4]
  960e57:	add    rax,0x1
  960e5b:	vaddss xmm15,xmm15,xmm14
  960e60:	vmovss xmm14,DWORD PTR [rsp+0x20]
  960e66:	vsubss xmm3,xmm3,xmm13
  960e6b:	vmulss xmm13,xmm7,xmm14
  960e70:	vmovss DWORD PTR [rsp+0x18],xmm15
  960e76:	vmulss xmm14,xmm6,xmm14
  960e7b:	vmulss xmm15,xmm7,xmm2
  960e7f:	vmulss xmm2,xmm6,xmm2
  960e83:	vaddss xmm15,xmm15,xmm14
  960e88:	vmovss xmm14,DWORD PTR [rsp+0x1c]
  960e8e:	vaddss xmm2,xmm2,xmm12
  960e93:	vmovss DWORD PTR [rsp+0x20],xmm15
  960e99:	vmulss xmm15,xmm1,xmm10
  960e9e:	vsubss xmm2,xmm2,xmm13
  960ea3:	vmulss xmm1,xmm1,xmm9
  960ea8:	vmulss xmm13,xmm10,xmm14
  960ead:	vmulss xmm14,xmm9,xmm14
  960eb2:	vaddss xmm1,xmm1,xmm12
  960eb7:	vaddss xmm15,xmm15,xmm14
  960ebc:	vmovss xmm14,DWORD PTR [rsp+0x10]
  960ec2:	vsubss xmm1,xmm1,xmm13
  960ec7:	vmulss xmm13,xmm8,xmm14
  960ecc:	vmovss DWORD PTR [rsp+0x1c],xmm15
  960ed2:	vmulss xmm14,xmm11,xmm14
  960ed7:	vmulss xmm15,xmm0,xmm8
  960edc:	vmulss xmm0,xmm0,xmm11
  960ee1:	vaddss xmm15,xmm15,xmm14
  960ee6:	vaddss xmm0,xmm0,xmm12
  960eeb:	vmovss DWORD PTR [rsp+0x10],xmm15
  960ef1:	vsubss xmm0,xmm0,xmm13
  960ef6:	cmp    rbx,rax
  960ef9:	jne    960dd8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x1b8>
  960eff:	vmovss xmm6,DWORD PTR [rsp+0x24]
  960f05:	vmovss DWORD PTR [r9+rdx*4],xmm5
  960f0b:	vmovss xmm7,DWORD PTR [rsp]
  960f10:	vmovss DWORD PTR [r9+rcx*1],xmm6
  960f16:	vmovss xmm6,DWORD PTR [rsp+0x18]
  960f1c:	vmovss DWORD PTR [r10+rdx*4],xmm4
  960f22:	mov    rax,QWORD PTR [rsp+0x28]
  960f27:	vmovss DWORD PTR [r10+rcx*1],xmm7
  960f2d:	vmovss xmm7,DWORD PTR [rsp+0x20]
  960f33:	vmovss DWORD PTR [r8+rdx*4],xmm3
  960f39:	vmovss DWORD PTR [r8+rcx*1],xmm6
  960f3f:	vmovss xmm6,DWORD PTR [rsp+0x1c]
  960f45:	vmovss DWORD PTR [r12+rdx*4],xmm2
  960f4b:	vmovss DWORD PTR [r12+rcx*1],xmm7
  960f51:	vmovss xmm7,DWORD PTR [rsp+0x10]
  960f57:	vmovss DWORD PTR [rbp+rdx*4+0x0],xmm1
  960f5d:	vmovss DWORD PTR [rbp+rcx*1+0x0],xmm6
  960f63:	vmovss DWORD PTR [r11+rdx*4],xmm0
  960f69:	add    rdx,0x1
  960f6d:	vmovss DWORD PTR [r11+rcx*1],xmm7
  960f73:	add    rcx,0x4
  960f77:	cmp    rdx,rax
  960f7a:	jne    960d38 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x118>
  960f80:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  960f8c:	lea    rbp,[rax-0x360]
  960f93:	mov    r12,rax
  960f96:	mov    r13,QWORD PTR [rbp+0x188]
  960f9d:	test   r13,r13
  960fa0:	je     960fc0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3a0>
  960fa2:	movzx  eax,BYTE PTR [rbp+0x180]
  960fa9:	lea    ebx,[rax-0x1]
  960fac:	mov    BYTE PTR [rbp+0x180],bl
  960fb2:	cmp    bl,0xf
  960fb5:	jbe    960fcf <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3af>
  960fb7:	add    WORD PTR [r13+0xfe2],0x1
  960fc0:	add    rsp,0x58
  960fc4:	pop    rbx
  960fc5:	pop    rbp
  960fc6:	pop    r12
  960fc8:	pop    r13
  960fca:	pop    r14
  960fcc:	pop    r15
  960fce:	ret
  960fcf:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  960fd4:	movzx  r15d,bl
  960fd8:	lea    rdx,[r15+r15*2]
  960fdc:	mov    rcx,rax
  960fdf:	vmovq  xmm6,rax
  960fe4:	shl    rdx,0x3
  960fe8:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  960fed:	vpinsrq xmm0,xmm6,rcx,0x1
  960ff3:	vpsubq xmm0,xmm0,XMMWORD PTR [r12+rdx*1-0x360]
  960ffd:	movzx  r12d,BYTE PTR [rbp+0x181]
  961005:	cmp    r12b,bl
  961008:	jae    96101e <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3fe>
  96100a:	lea    eax,[r15-0x1]
  96100e:	vmovq  rdx,xmm0
  961013:	cdqe
  961015:	lea    rax,[rax+rax*2]
  961019:	add    QWORD PTR [rbp+rax*8+0x8],rdx
  96101e:	movzx  ebp,WORD PTR [r13+0xfe0]
  961026:	cmp    bp,0x7f
  96102a:	je     96117b <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x55b>
  961030:	vmovdqa XMMWORD PTR [rsp],xmm0
  961035:	sub    ebx,r12d
  961038:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  961044:	lea    rdx,[r15+r15*2]
  961048:	vmovdqa xmm0,XMMWORD PTR [rsp]
  96104d:	mov    DWORD PTR [rsp+0x48],0x0
  961055:	vmovdqu XMMWORD PTR [rsp+0x38],xmm0
  96105b:	lea    rdx,[rax+rdx*8-0x360]
  961063:	mov    rax,QWORD PTR [rdx]
  961066:	mov    QWORD PTR [rsp+0x30],rax
  96106b:	movzx  eax,BYTE PTR [rdx+0x10]
  96106f:	lea    edx,[rbp+0x1]
  961072:	vmovdqa xmm7,XMMWORD PTR [rsp+0x30]
  961078:	mov    WORD PTR [r13+0xfe0],dx
  961080:	mov    ah,bl
  961082:	mov    WORD PTR [rsp+0x4c],ax
  961087:	movzx  eax,bp
  96108a:	vmovdqu xmm6,XMMWORD PTR [rsp+0x3e]
  961090:	shl    rax,0x5
  961094:	add    rax,r13
  961097:	vmovdqu XMMWORD PTR [rax],xmm7
  96109b:	vmovdqu XMMWORD PTR [rax+0xe],xmm6
  9610a0:	add    rsp,0x58
  9610a4:	pop    rbx
  9610a5:	pop    rbp
  9610a6:	pop    r12
  9610a8:	pop    r13
  9610aa:	pop    r14
  9610ac:	pop    r15
  9610ae:	ret
  9610af:	shl    rdi,0x4
  9610b3:	add    rsi,rdi
  9610b6:	cmp    r9,rsi
  9610b9:	je     960cd4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb4>
  9610bf:	mov    QWORD PTR [r12+0x8],rsi
  9610c4:	jmp    960cd4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb4>
  9610c9:	mov    eax,0x1
  9610ce:	lock xadd DWORD PTR [rip+0x2c5a8a],eax        # c26b60 <cpl::Profiling::registerRegion(char const*)::counter>
  9610d6:	mov    ecx,0x1
  9610db:	mov    ebp,0x1
  9610e0:	add    eax,0x2
  9610e3:	cmp    eax,0xfe
  9610e8:	ja     961102 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4e2>
  9610ea:	lea    rdx,[rip+0x2c3f2f]        # c25020 <cpl::Profiling::regions>
  9610f1:	mov    ecx,eax
  9610f3:	mov    ebp,eax
  9610f5:	lea    rdi,[rip+0x9f1b6]        # a002b2 <_IO_stdin_used+0x182b2>
  9610fc:	mov    QWORD PTR [rdx+rcx*8],rdi
  961100:	mov    ecx,eax
  961102:	xor    eax,eax
  961104:	lock cmpxchg BYTE PTR [rip+0x2d0fd1],cl        # c320dd <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  96110c:	cmovne ebp,eax
  96110f:	jmp    960c4f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x2f>
  961114:	mov    rsi,rcx
  961117:	mov    rdi,r12
  96111a:	sub    rsi,rax
  96111d:	call   939b50 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  961122:	mov    rdx,QWORD PTR [r13+0x38]
  961126:	mov    r8,QWORD PTR [r13+0x48]
  96112a:	jmp    960cd4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb4>
  96112f:	mov    QWORD PTR [rsp+0x10],rcx
  961134:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  961139:	mov    rcx,QWORD PTR [rsp+0x10]
  96113e:	mov    BYTE PTR [rsp+0x40],bpl
  961143:	mov    QWORD PTR [rsp+0x30],rax
  961148:	lea    rax,[r15+r15*2]
  96114c:	mov    QWORD PTR [rsp+0x38],0x0
  961155:	vmovdqa xmm7,XMMWORD PTR [rsp+0x30]
  96115b:	shl    rax,0x3
  96115f:	vmovdqu XMMWORD PTR [r14+rax*1-0x360],xmm7
  961169:	mov    BYTE PTR [rax+rcx*1+0x10],bpl
  96116e:	movzx  r15d,BYTE PTR [rcx+0x180]
  961176:	jmp    960c81 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x61>
  96117b:	movzx  eax,WORD PTR [r13+0xfe2]
  961183:	cmp    ax,0xffff
  961187:	je     960fc0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3a0>
  96118d:	add    eax,0x1
  961190:	mov    WORD PTR [r13+0xfe2],ax
  961198:	jmp    960fc0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3a0>
  96119d:	endbr64
  9611a1:	mov    rbx,rax
  9611a4:	vzeroupper
  9611a7:	call   93c730 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  9611ac:	mov    rdi,rbx
  9611af:	call   e1b10 <_Unwind_Resume@plt>

Disassembly of section .fini:
