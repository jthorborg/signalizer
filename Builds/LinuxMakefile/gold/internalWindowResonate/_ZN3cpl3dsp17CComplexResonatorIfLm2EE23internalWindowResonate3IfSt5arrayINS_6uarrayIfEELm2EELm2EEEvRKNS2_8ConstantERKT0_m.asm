; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000965b00 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long):
  965b00:	endbr64
  965b04:	push   r15
  965b06:	push   r14
  965b08:	push   r13
  965b0a:	mov    r13,rsi
  965b0d:	push   r12
  965b0f:	mov    r12,rdi
  965b12:	push   rbp
  965b13:	push   rbx
  965b14:	mov    rbx,rcx
  965b17:	sub    rsp,0x58
  965b1b:	mov    QWORD PTR [rsp],rdx
  965b1f:	movzx  ebp,BYTE PTR [rip+0x2cd5ad]        # c330d3 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  965b26:	test   bpl,bpl
  965b29:	je     965fa9 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4a9>
  965b2f:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  965b3b:	cmp    QWORD PTR [rax-0x1d8],0x0
  965b43:	mov    r14,rax
  965b46:	lea    rcx,[rax-0x360]
  965b4d:	je     965b78 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x78>
  965b4f:	movzx  r15d,BYTE PTR [rax-0x1e0]
  965b57:	cmp    r15b,0xf
  965b5b:	jbe    96600f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x50f>
  965b61:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  965b6d:	add    r15d,0x1
  965b71:	mov    BYTE PTR [rax-0x1e0],r15b
  965b78:	mov    rdx,QWORD PTR [r13+0x38]
  965b7c:	mov    r8,QWORD PTR [r13+0x48]
  965b80:	mov    r9,QWORD PTR [r12+0x8]
  965b85:	mov    rsi,QWORD PTR [r12]
  965b89:	mov    rdi,rdx
  965b8c:	imul   rdi,r8
  965b90:	mov    rax,r9
  965b93:	sub    rax,rsi
  965b96:	sar    rax,0x2
  965b9a:	lea    rcx,[rdi*4+0x0]
  965ba2:	cmp    rax,rcx
  965ba5:	jb     965ff4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4f4>
  965bab:	cmp    rcx,rax
  965bae:	jb     965f8f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x48f>
  965bb4:	mov    rsi,QWORD PTR [r13+0x40]
  965bb8:	lea    rax,[r8+r8*1]
  965bbc:	imul   rdx,rax
  965bc0:	mov    QWORD PTR [rsp+0x28],rsi
  965bc5:	test   rsi,rsi
  965bc8:	je     965e60 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x360>
  965bce:	mov    rsi,QWORD PTR [rsp]
  965bd2:	lea    rcx,[r8*4+0x0]
  965bda:	shl    r8,0x3
  965bde:	mov    r9,QWORD PTR [r12]
  965be2:	lea    r11,[rcx+rdx*1]
  965be6:	mov    r14,QWORD PTR [r13+0x0]
  965bea:	mov    rdi,QWORD PTR [rsi]
  965bed:	mov    rsi,QWORD PTR [rsi+0x10]
  965bf1:	lea    r12,[r9+rdx*4]
  965bf5:	mov    rdx,r11
  965bf8:	lea    r13,[r14+r8*1]
  965bfc:	sub    rdx,rax
  965bff:	lea    r10,[r9+r8*1]
  965c03:	lea    r15,[r13+r8*1+0x0]
  965c08:	add    r8,r10
  965c0b:	lea    rbp,[r9+rdx*4]
  965c0f:	xor    edx,edx
  965c11:	lea    r11,[r9+r11*4]
  965c15:	nop    DWORD PTR [rax]
  965c18:	vmovss xmm1,DWORD PTR [r9+rcx*1]
  965c1e:	vmovss xmm0,DWORD PTR [r8+rcx*1]
  965c24:	vmovss xmm2,DWORD PTR [r10+rcx*1]
  965c2a:	vmovss DWORD PTR [rsp+0x24],xmm1
  965c30:	vmovss xmm12,DWORD PTR [r11+rcx*1]
  965c36:	vmovss DWORD PTR [rsp+0x18],xmm0
  965c3c:	vmovss xmm1,DWORD PTR [r12+rcx*1]
  965c42:	vmovss DWORD PTR [rsp],xmm2
  965c47:	vmovss xmm0,DWORD PTR [rbp+rcx*1+0x0]
  965c4d:	vmovss DWORD PTR [rsp+0x10],xmm12
  965c53:	vmovss DWORD PTR [rsp+0x20],xmm1
  965c59:	vmovss xmm6,DWORD PTR [r14+rdx*4]
  965c5f:	vmovss DWORD PTR [rsp+0x1c],xmm0
  965c65:	vmovss xmm7,DWORD PTR [r14+rcx*1]
  965c6b:	vmovss xmm9,DWORD PTR [r13+rdx*4+0x0]
  965c72:	vmovss xmm10,DWORD PTR [r13+rcx*1+0x0]
  965c79:	vmovss xmm11,DWORD PTR [r15+rdx*4]
  965c7f:	vmovss xmm8,DWORD PTR [r15+rcx*1]
  965c85:	vmovss xmm5,DWORD PTR [r9+rdx*4]
  965c8b:	vmovss xmm4,DWORD PTR [r10+rdx*4]
  965c91:	vmovss xmm3,DWORD PTR [r8+rdx*4]
  965c97:	vmovss xmm2,DWORD PTR [r12+rdx*4]
  965c9d:	vmovss xmm1,DWORD PTR [rbp+rdx*4+0x0]
  965ca3:	vmovss xmm0,DWORD PTR [r11+rdx*4]
  965ca9:	test   rbx,rbx
  965cac:	je     965ddf <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x2df>
  965cb2:	xor    eax,eax
  965cb4:	nop    DWORD PTR [rax+0x0]
  965cb8:	vmovss xmm15,DWORD PTR [rsp+0x24]
  965cbe:	vmulss xmm14,xmm5,xmm7
  965cc2:	vmulss xmm5,xmm5,xmm6
  965cc6:	vmovss xmm12,DWORD PTR [rdi+rax*4]
  965ccb:	vmulss xmm13,xmm7,xmm15
  965cd0:	vmulss xmm15,xmm6,xmm15
  965cd5:	vaddss xmm5,xmm5,xmm12
  965cda:	vaddss xmm15,xmm14,xmm15
  965cdf:	vmovss xmm14,DWORD PTR [rsp]
  965ce4:	vsubss xmm5,xmm5,xmm13
  965ce9:	vmulss xmm13,xmm10,xmm14
  965cee:	vmovss DWORD PTR [rsp+0x24],xmm15
  965cf4:	vmulss xmm14,xmm9,xmm14
  965cf9:	vmulss xmm15,xmm10,xmm4
  965cfd:	vmulss xmm4,xmm9,xmm4
  965d01:	vaddss xmm15,xmm15,xmm14
  965d06:	vmovss xmm14,DWORD PTR [rsp+0x18]
  965d0c:	vaddss xmm4,xmm4,xmm12
  965d11:	vmovss DWORD PTR [rsp],xmm15
  965d16:	vmulss xmm15,xmm8,xmm3
  965d1a:	vsubss xmm4,xmm4,xmm13
  965d1f:	vmulss xmm3,xmm11,xmm3
  965d23:	vmulss xmm13,xmm8,xmm14
  965d28:	vmulss xmm14,xmm11,xmm14
  965d2d:	vaddss xmm3,xmm3,xmm12
  965d32:	vmovss xmm12,DWORD PTR [rsi+rax*4]
  965d37:	add    rax,0x1
  965d3b:	vaddss xmm15,xmm15,xmm14
  965d40:	vmovss xmm14,DWORD PTR [rsp+0x20]
  965d46:	vsubss xmm3,xmm3,xmm13
  965d4b:	vmulss xmm13,xmm7,xmm14
  965d50:	vmovss DWORD PTR [rsp+0x18],xmm15
  965d56:	vmulss xmm14,xmm6,xmm14
  965d5b:	vmulss xmm15,xmm7,xmm2
  965d5f:	vmulss xmm2,xmm6,xmm2
  965d63:	vaddss xmm15,xmm15,xmm14
  965d68:	vmovss xmm14,DWORD PTR [rsp+0x1c]
  965d6e:	vaddss xmm2,xmm2,xmm12
  965d73:	vmovss DWORD PTR [rsp+0x20],xmm15
  965d79:	vmulss xmm15,xmm1,xmm10
  965d7e:	vsubss xmm2,xmm2,xmm13
  965d83:	vmulss xmm1,xmm1,xmm9
  965d88:	vmulss xmm13,xmm10,xmm14
  965d8d:	vmulss xmm14,xmm9,xmm14
  965d92:	vaddss xmm1,xmm1,xmm12
  965d97:	vaddss xmm15,xmm15,xmm14
  965d9c:	vmovss xmm14,DWORD PTR [rsp+0x10]
  965da2:	vsubss xmm1,xmm1,xmm13
  965da7:	vmulss xmm13,xmm8,xmm14
  965dac:	vmovss DWORD PTR [rsp+0x1c],xmm15
  965db2:	vmulss xmm14,xmm11,xmm14
  965db7:	vmulss xmm15,xmm0,xmm8
  965dbc:	vmulss xmm0,xmm0,xmm11
  965dc1:	vaddss xmm15,xmm15,xmm14
  965dc6:	vaddss xmm0,xmm0,xmm12
  965dcb:	vmovss DWORD PTR [rsp+0x10],xmm15
  965dd1:	vsubss xmm0,xmm0,xmm13
  965dd6:	cmp    rbx,rax
  965dd9:	jne    965cb8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x1b8>
  965ddf:	vmovss xmm6,DWORD PTR [rsp+0x24]
  965de5:	vmovss DWORD PTR [r9+rdx*4],xmm5
  965deb:	vmovss xmm7,DWORD PTR [rsp]
  965df0:	vmovss DWORD PTR [r9+rcx*1],xmm6
  965df6:	vmovss xmm6,DWORD PTR [rsp+0x18]
  965dfc:	vmovss DWORD PTR [r10+rdx*4],xmm4
  965e02:	mov    rax,QWORD PTR [rsp+0x28]
  965e07:	vmovss DWORD PTR [r10+rcx*1],xmm7
  965e0d:	vmovss xmm7,DWORD PTR [rsp+0x20]
  965e13:	vmovss DWORD PTR [r8+rdx*4],xmm3
  965e19:	vmovss DWORD PTR [r8+rcx*1],xmm6
  965e1f:	vmovss xmm6,DWORD PTR [rsp+0x1c]
  965e25:	vmovss DWORD PTR [r12+rdx*4],xmm2
  965e2b:	vmovss DWORD PTR [r12+rcx*1],xmm7
  965e31:	vmovss xmm7,DWORD PTR [rsp+0x10]
  965e37:	vmovss DWORD PTR [rbp+rdx*4+0x0],xmm1
  965e3d:	vmovss DWORD PTR [rbp+rcx*1+0x0],xmm6
  965e43:	vmovss DWORD PTR [r11+rdx*4],xmm0
  965e49:	add    rdx,0x1
  965e4d:	vmovss DWORD PTR [r11+rcx*1],xmm7
  965e53:	add    rcx,0x4
  965e57:	cmp    rdx,rax
  965e5a:	jne    965c18 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x118>
  965e60:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  965e6c:	lea    rbp,[rax-0x360]
  965e73:	mov    r12,rax
  965e76:	mov    r13,QWORD PTR [rbp+0x188]
  965e7d:	test   r13,r13
  965e80:	je     965ea0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3a0>
  965e82:	movzx  eax,BYTE PTR [rbp+0x180]
  965e89:	lea    ebx,[rax-0x1]
  965e8c:	mov    BYTE PTR [rbp+0x180],bl
  965e92:	cmp    bl,0xf
  965e95:	jbe    965eaf <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3af>
  965e97:	add    WORD PTR [r13+0xfe2],0x1
  965ea0:	add    rsp,0x58
  965ea4:	pop    rbx
  965ea5:	pop    rbp
  965ea6:	pop    r12
  965ea8:	pop    r13
  965eaa:	pop    r14
  965eac:	pop    r15
  965eae:	ret
  965eaf:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  965eb4:	movzx  r15d,bl
  965eb8:	lea    rdx,[r15+r15*2]
  965ebc:	mov    rcx,rax
  965ebf:	vmovq  xmm6,rax
  965ec4:	shl    rdx,0x3
  965ec8:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  965ecd:	vpinsrq xmm0,xmm6,rcx,0x1
  965ed3:	vpsubq xmm0,xmm0,XMMWORD PTR [r12+rdx*1-0x360]
  965edd:	movzx  r12d,BYTE PTR [rbp+0x181]
  965ee5:	cmp    r12b,bl
  965ee8:	jae    965efe <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3fe>
  965eea:	lea    eax,[r15-0x1]
  965eee:	vmovq  rdx,xmm0
  965ef3:	cdqe
  965ef5:	lea    rax,[rax+rax*2]
  965ef9:	add    QWORD PTR [rbp+rax*8+0x8],rdx
  965efe:	movzx  ebp,WORD PTR [r13+0xfe0]
  965f06:	cmp    bp,0x7f
  965f0a:	je     96605b <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x55b>
  965f10:	vmovdqa XMMWORD PTR [rsp],xmm0
  965f15:	sub    ebx,r12d
  965f18:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  965f24:	lea    rdx,[r15+r15*2]
  965f28:	vmovdqa xmm0,XMMWORD PTR [rsp]
  965f2d:	mov    DWORD PTR [rsp+0x48],0x0
  965f35:	vmovdqu XMMWORD PTR [rsp+0x38],xmm0
  965f3b:	lea    rdx,[rax+rdx*8-0x360]
  965f43:	mov    rax,QWORD PTR [rdx]
  965f46:	mov    QWORD PTR [rsp+0x30],rax
  965f4b:	movzx  eax,BYTE PTR [rdx+0x10]
  965f4f:	lea    edx,[rbp+0x1]
  965f52:	vmovdqa xmm7,XMMWORD PTR [rsp+0x30]
  965f58:	mov    WORD PTR [r13+0xfe0],dx
  965f60:	mov    ah,bl
  965f62:	mov    WORD PTR [rsp+0x4c],ax
  965f67:	movzx  eax,bp
  965f6a:	vmovdqu xmm6,XMMWORD PTR [rsp+0x3e]
  965f70:	shl    rax,0x5
  965f74:	add    rax,r13
  965f77:	vmovdqu XMMWORD PTR [rax],xmm7
  965f7b:	vmovdqu XMMWORD PTR [rax+0xe],xmm6
  965f80:	add    rsp,0x58
  965f84:	pop    rbx
  965f85:	pop    rbp
  965f86:	pop    r12
  965f88:	pop    r13
  965f8a:	pop    r14
  965f8c:	pop    r15
  965f8e:	ret
  965f8f:	shl    rdi,0x4
  965f93:	add    rsi,rdi
  965f96:	cmp    r9,rsi
  965f99:	je     965bb4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb4>
  965f9f:	mov    QWORD PTR [r12+0x8],rsi
  965fa4:	jmp    965bb4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb4>
  965fa9:	mov    eax,0x1
  965fae:	lock xadd DWORD PTR [rip+0x2c1baa],eax        # c27b60 <cpl::Profiling::registerRegion(char const*)::counter>
  965fb6:	mov    ecx,0x1
  965fbb:	mov    ebp,0x1
  965fc0:	add    eax,0x2
  965fc3:	cmp    eax,0xfe
  965fc8:	ja     965fe2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4e2>
  965fca:	lea    rdx,[rip+0x2c004f]        # c26020 <cpl::Profiling::regions>
  965fd1:	mov    ecx,eax
  965fd3:	mov    ebp,eax
  965fd5:	lea    rdi,[rip+0x9b2d6]        # a012b2 <_IO_stdin_used+0x182b2>
  965fdc:	mov    QWORD PTR [rdx+rcx*8],rdi
  965fe0:	mov    ecx,eax
  965fe2:	xor    eax,eax
  965fe4:	lock cmpxchg BYTE PTR [rip+0x2cd0e7],cl        # c330d3 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  965fec:	cmovne ebp,eax
  965fef:	jmp    965b2f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x2f>
  965ff4:	mov    rsi,rcx
  965ff7:	mov    rdi,r12
  965ffa:	sub    rsi,rax
  965ffd:	call   939bb0 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  966002:	mov    rdx,QWORD PTR [r13+0x38]
  966006:	mov    r8,QWORD PTR [r13+0x48]
  96600a:	jmp    965bb4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb4>
  96600f:	mov    QWORD PTR [rsp+0x10],rcx
  966014:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  966019:	mov    rcx,QWORD PTR [rsp+0x10]
  96601e:	mov    BYTE PTR [rsp+0x40],bpl
  966023:	mov    QWORD PTR [rsp+0x30],rax
  966028:	lea    rax,[r15+r15*2]
  96602c:	mov    QWORD PTR [rsp+0x38],0x0
  966035:	vmovdqa xmm7,XMMWORD PTR [rsp+0x30]
  96603b:	shl    rax,0x3
  96603f:	vmovdqu XMMWORD PTR [r14+rax*1-0x360],xmm7
  966049:	mov    BYTE PTR [rax+rcx*1+0x10],bpl
  96604e:	movzx  r15d,BYTE PTR [rcx+0x180]
  966056:	jmp    965b61 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x61>
  96605b:	movzx  eax,WORD PTR [r13+0xfe2]
  966063:	cmp    ax,0xffff
  966067:	je     965ea0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3a0>
  96606d:	add    eax,0x1
  966070:	mov    WORD PTR [r13+0xfe2],ax
  966078:	jmp    965ea0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3a0>
  96607d:	endbr64
  966081:	mov    rbx,rax
  966084:	vzeroupper
  966087:	call   93c790 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  96608c:	mov    rdi,rbx
  96608f:	call   e1b10 <_Unwind_Resume@plt>

Disassembly of section .fini:
