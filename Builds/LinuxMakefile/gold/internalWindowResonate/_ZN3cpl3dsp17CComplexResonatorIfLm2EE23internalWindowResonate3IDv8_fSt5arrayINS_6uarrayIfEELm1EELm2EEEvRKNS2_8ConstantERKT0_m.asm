; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

000000000094eaa0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long):
  94eaa0:	endbr64
  94eaa4:	push   rbp
  94eaa5:	mov    rbp,rsp
  94eaa8:	push   r15
  94eaaa:	push   r14
  94eaac:	mov    r14,rdi
  94eaaf:	push   r13
  94eab1:	push   r12
  94eab3:	push   rbx
  94eab4:	mov    rbx,rcx
  94eab7:	and    rsp,0xffffffffffffffe0
  94eabb:	sub    rsp,0x180
  94eac2:	mov    QWORD PTR [rsp+0x90],rsi
  94eaca:	mov    QWORD PTR [rsp+0x40],rdx
  94eacf:	movzx  r13d,BYTE PTR [rip+0x2e4631]        # c33108 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  94ead7:	test   r13b,r13b
  94eada:	je     94f171 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x6d1>
  94eae0:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94eaec:	cmp    QWORD PTR [rax-0x1d8],0x0
  94eaf4:	lea    r15,[rax-0x360]
  94eafb:	je     94eb26 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x86>
  94eafd:	movzx  r12d,BYTE PTR [r15+0x180]
  94eb05:	cmp    r12b,0xf
  94eb09:	jbe    94f1df <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x73f>
  94eb0f:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94eb1b:	add    r12d,0x1
  94eb1f:	mov    BYTE PTR [rax-0x1e0],r12b
  94eb26:	mov    rax,QWORD PTR [rsp+0x90]
  94eb2e:	mov    r9,QWORD PTR [r14+0x8]
  94eb32:	mov    rdi,QWORD PTR [r14]
  94eb35:	mov    rdx,QWORD PTR [rax+0x48]
  94eb39:	mov    rax,QWORD PTR [rax+0x38]
  94eb3d:	mov    rcx,r9
  94eb40:	sub    rcx,rdi
  94eb43:	sar    rcx,0x2
  94eb47:	mov    r8,rax
  94eb4a:	imul   r8,rdx
  94eb4e:	lea    rsi,[r8*4+0x0]
  94eb56:	cmp    rcx,rsi
  94eb59:	jb     94f1bf <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x71f>
  94eb5f:	cmp    rsi,rcx
  94eb62:	jb     94f158 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x6b8>
  94eb68:	mov    rdi,QWORD PTR [rsp+0x90]
  94eb70:	lea    rcx,[rdx+rdx*1]
  94eb74:	imul   rax,rcx
  94eb78:	cmp    QWORD PTR [rdi+0x40],0x0
  94eb7d:	je     94eff3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x553>
  94eb83:	lea    r11,[rdx*4+0x0]
  94eb8b:	xor    r10d,r10d
  94eb8e:	mov    QWORD PTR [rsp+0x98],0x0
  94eb9a:	mov    r12,rbx
  94eb9d:	lea    rsi,[r11+rax*1]
  94eba1:	shl    rax,0x2
  94eba5:	mov    QWORD PTR [rsp+0x38],r11
  94ebaa:	mov    QWORD PTR [rsp+0x20],rax
  94ebaf:	lea    rdi,[rdx*8+0x0]
  94ebb7:	mov    rax,rsi
  94ebba:	shl    rdx,0x4
  94ebbe:	sub    rax,rcx
  94ebc1:	mov    QWORD PTR [rsp+0x30],rdi
  94ebc6:	shl    rax,0x2
  94ebca:	mov    QWORD PTR [rsp+0x28],rdx
  94ebcf:	mov    QWORD PTR [rsp+0x18],rax
  94ebd4:	lea    rax,[rsi*4+0x0]
  94ebdc:	mov    QWORD PTR [rsp+0x88],rax
  94ebe4:	nop    DWORD PTR [rax+0x0]
  94ebe8:	mov    rax,QWORD PTR [rsp+0x90]
  94ebf0:	mov    rdi,QWORD PTR [rsp+0x30]
  94ebf5:	mov    rbx,QWORD PTR [rsp+0x38]
  94ebfa:	mov    rcx,QWORD PTR [rsp+0x28]
  94ebff:	mov    rsi,QWORD PTR [r14]
  94ec02:	mov    rax,QWORD PTR [rax]
  94ec05:	lea    r8,[r10+rdi*1]
  94ec09:	mov    rdx,QWORD PTR [rsp+0x40]
  94ec0e:	lea    rdi,[rbx+r8*1]
  94ec12:	lea    r9,[r11+rcx*1]
  94ec16:	add    rbx,rdi
  94ec19:	vmovups xmm1,XMMWORD PTR [rax+r8*1]
  94ec1f:	vinsertf128 ymm9,ymm1,XMMWORD PTR [rax+r8*1+0x10],0x1
  94ec27:	vmovups xmm2,XMMWORD PTR [rax+rdi*1]
  94ec2c:	vinsertf128 ymm10,ymm2,XMMWORD PTR [rax+rdi*1+0x10],0x1
  94ec34:	vmovups xmm5,XMMWORD PTR [rax+rbx*1]
  94ec39:	vinsertf128 ymm11,ymm5,XMMWORD PTR [rax+rbx*1+0x10],0x1
  94ec41:	vmovups xmm1,XMMWORD PTR [rax+r9*1]
  94ec47:	vinsertf128 ymm8,ymm1,XMMWORD PTR [rax+r9*1+0x10],0x1
  94ec4f:	vmovups xmm6,XMMWORD PTR [rax+r10*1]
  94ec55:	vinsertf128 ymm6,ymm6,XMMWORD PTR [rax+r10*1+0x10],0x1
  94ec5d:	vmovups xmm7,XMMWORD PTR [rax+r11*1]
  94ec63:	vinsertf128 ymm7,ymm7,XMMWORD PTR [rax+r11*1+0x10],0x1
  94ec6b:	lea    rax,[rsi+r10*1]
  94ec6f:	vmovups xmm2,XMMWORD PTR [rax]
  94ec73:	vinsertf128 ymm0,ymm2,XMMWORD PTR [rax+0x10],0x1
  94ec7a:	mov    rcx,QWORD PTR [rdx]
  94ec7d:	vmovups xmm2,XMMWORD PTR [rsi+r8*1]
  94ec83:	vinsertf128 ymm1,ymm2,XMMWORD PTR [rsi+r8*1+0x10],0x1
  94ec8b:	vmovaps YMMWORD PTR [rsp+0x120],ymm1
  94ec94:	vmovups xmm2,XMMWORD PTR [rsi+rdi*1]
  94ec99:	vinsertf128 ymm1,ymm2,XMMWORD PTR [rsi+rdi*1+0x10],0x1
  94eca1:	vmovaps YMMWORD PTR [rsp+0x100],ymm1
  94ecaa:	vmovups xmm2,XMMWORD PTR [rsi+rbx*1]
  94ecaf:	vinsertf128 ymm1,ymm2,XMMWORD PTR [rsi+rbx*1+0x10],0x1
  94ecb7:	vmovups xmm5,XMMWORD PTR [rsi+r11*1]
  94ecbd:	vinsertf128 ymm5,ymm5,XMMWORD PTR [rsi+r11*1+0x10],0x1
  94ecc5:	vmovaps YMMWORD PTR [rsp+0x140],ymm1
  94ecce:	mov    rdx,QWORD PTR [rdx+0x10]
  94ecd2:	vmovups xmm2,XMMWORD PTR [rsi+r9*1]
  94ecd8:	vinsertf128 ymm1,ymm2,XMMWORD PTR [rsi+r9*1+0x10],0x1
  94ece0:	vmovaps YMMWORD PTR [rsp+0xc0],ymm1
  94ece9:	mov    QWORD PTR [rsp+0x48],rdx
  94ecee:	mov    rdx,QWORD PTR [rsp+0x20]
  94ecf3:	lea    r15,[r10+rdx*1]
  94ecf7:	vmovups xmm1,XMMWORD PTR [rsi+r15*1]
  94ecfd:	mov    QWORD PTR [rsp+0x80],r15
  94ed05:	vinsertf128 ymm4,ymm1,XMMWORD PTR [rsi+r15*1+0x10],0x1
  94ed0d:	lea    r15,[r11+rdx*1]
  94ed11:	mov    rdx,QWORD PTR [rsp+0x18]
  94ed16:	mov    QWORD PTR [rsp+0x78],r15
  94ed1b:	vmovups xmm1,XMMWORD PTR [rsi+r15*1]
  94ed21:	vinsertf128 ymm2,ymm1,XMMWORD PTR [rsi+r15*1+0x10],0x1
  94ed29:	vmovaps YMMWORD PTR [rsp+0xa0],ymm2
  94ed32:	lea    r15,[r10+rdx*1]
  94ed36:	lea    r13,[r11+rdx*1]
  94ed3a:	mov    rdx,QWORD PTR [rsp+0x88]
  94ed42:	vmovups xmm2,XMMWORD PTR [rsi+r15*1]
  94ed48:	vinsertf128 ymm3,ymm2,XMMWORD PTR [rsi+r15*1+0x10],0x1
  94ed50:	vmovups xmm2,XMMWORD PTR [rsi+r13*1]
  94ed56:	vinsertf128 ymm1,ymm2,XMMWORD PTR [rsi+r13*1+0x10],0x1
  94ed5e:	vmovaps YMMWORD PTR [rsp+0xe0],ymm1
  94ed67:	add    rdx,r10
  94ed6a:	vmovups xmm1,XMMWORD PTR [rsi+rdx*1]
  94ed6f:	mov    QWORD PTR [rsp+0x70],rdx
  94ed74:	vinsertf128 ymm2,ymm1,XMMWORD PTR [rsi+rdx*1+0x10],0x1
  94ed7c:	mov    rdx,QWORD PTR [rsp+0x88]
  94ed84:	vmovaps XMMWORD PTR [rsp+0x60],xmm1
  94ed8a:	add    rdx,r11
  94ed8d:	vmovups xmm1,XMMWORD PTR [rsi+rdx*1]
  94ed92:	mov    QWORD PTR [rsp+0x60],rdx
  94ed97:	vmovaps XMMWORD PTR [rsp+0x50],xmm1
  94ed9d:	vinsertf128 ymm1,ymm1,XMMWORD PTR [rsi+rdx*1+0x10],0x1
  94eda5:	test   r12,r12
  94eda8:	je     94ef17 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x477>
  94edae:	mov    rdx,QWORD PTR [rsp+0x48]
  94edb3:	mov    rax,rcx
  94edb6:	xor    ecx,ecx
  94edb8:	nop    DWORD PTR [rax+rax*1+0x0]
  94edc0:	vmulps ymm14,ymm0,ymm7
  94edc4:	vbroadcastss ymm12,DWORD PTR [rax]
  94edc9:	vmulps ymm0,ymm0,ymm6
  94edcd:	vmulps ymm13,ymm7,ymm5
  94edd1:	vmulps ymm5,ymm6,ymm5
  94edd5:	vaddps ymm5,ymm14,ymm5
  94edd9:	add    rcx,0x1
  94eddd:	add    rdx,0x4
  94ede1:	vaddps ymm0,ymm0,ymm12
  94ede6:	add    rax,0x4
  94edea:	vmovaps ymm14,YMMWORD PTR [rsp+0x100]
  94edf3:	vmulps ymm15,ymm10,YMMWORD PTR [rsp+0x120]
  94edfc:	vsubps ymm0,ymm0,ymm13
  94ee01:	vmulps ymm13,ymm10,ymm14
  94ee06:	vmulps ymm14,ymm9,ymm14
  94ee0b:	vaddps ymm14,ymm15,ymm14
  94ee10:	vmulps ymm15,ymm8,YMMWORD PTR [rsp+0x140]
  94ee19:	vmovaps YMMWORD PTR [rsp+0x100],ymm14
  94ee22:	vmulps ymm14,ymm9,YMMWORD PTR [rsp+0x120]
  94ee2b:	vaddps ymm14,ymm14,ymm12
  94ee30:	vsubps ymm13,ymm14,ymm13
  94ee35:	vmovaps ymm14,YMMWORD PTR [rsp+0xc0]
  94ee3e:	vmovaps YMMWORD PTR [rsp+0x120],ymm13
  94ee47:	vmulps ymm13,ymm8,ymm14
  94ee4c:	vmulps ymm14,ymm11,ymm14
  94ee51:	vaddps ymm14,ymm15,ymm14
  94ee56:	vmovaps ymm15,YMMWORD PTR [rsp+0xa0]
  94ee5f:	vmovaps YMMWORD PTR [rsp+0xc0],ymm14
  94ee68:	vmulps ymm14,ymm11,YMMWORD PTR [rsp+0x140]
  94ee71:	vaddps ymm14,ymm14,ymm12
  94ee76:	vmulps ymm12,ymm6,ymm4
  94ee7a:	vmulps ymm4,ymm7,ymm4
  94ee7e:	vsubps ymm13,ymm14,ymm13
  94ee83:	vmulps ymm14,ymm7,ymm15
  94ee88:	vsubps ymm12,ymm12,ymm14
  94ee8d:	vmulps ymm14,ymm6,ymm15
  94ee92:	vmovaps ymm15,YMMWORD PTR [rsp+0xe0]
  94ee9b:	vaddps ymm4,ymm4,ymm14
  94eea0:	vmovaps YMMWORD PTR [rsp+0x140],ymm13
  94eea9:	vbroadcastss ymm13,DWORD PTR [rdx-0x4]
  94eeaf:	vmulps ymm14,ymm10,ymm15
  94eeb4:	vmovaps YMMWORD PTR [rsp+0xa0],ymm4
  94eebd:	vaddps ymm4,ymm13,ymm12
  94eec2:	vmulps ymm12,ymm3,ymm9
  94eec7:	vmulps ymm3,ymm3,ymm10
  94eecc:	vsubps ymm12,ymm12,ymm14
  94eed1:	vmulps ymm14,ymm9,ymm15
  94eed6:	vaddps ymm3,ymm3,ymm14
  94eedb:	vmulps ymm14,ymm2,ymm11
  94eee0:	vmulps ymm2,ymm8,ymm2
  94eee4:	vmovaps YMMWORD PTR [rsp+0xe0],ymm3
  94eeed:	vaddps ymm3,ymm13,ymm12
  94eef2:	vmulps ymm12,ymm1,ymm8
  94eef7:	vmulps ymm1,ymm1,ymm11
  94eefc:	vsubps ymm12,ymm14,ymm12
  94ef01:	vaddps ymm1,ymm2,ymm1
  94ef05:	vaddps ymm2,ymm13,ymm12
  94ef0a:	cmp    r12,rcx
  94ef0d:	jne    94edc0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x320>
  94ef13:	lea    rax,[rsi+r10*1]
  94ef17:	vmovaps YMMWORD PTR [rax],ymm0
  94ef1b:	mov    rax,QWORD PTR [r14]
  94ef1e:	add    r10,0x20
  94ef22:	vmovaps ymm6,YMMWORD PTR [rsp+0x120]
  94ef2b:	vmovaps ymm7,YMMWORD PTR [rsp+0x100]
  94ef34:	vmovaps YMMWORD PTR [rax+r11*1],ymm5
  94ef3a:	mov    rax,QWORD PTR [r14]
  94ef3d:	add    r11,0x20
  94ef41:	vmovaps YMMWORD PTR [rax+r8*1],ymm6
  94ef47:	mov    rax,QWORD PTR [r14]
  94ef4a:	vmovaps ymm6,YMMWORD PTR [rsp+0x140]
  94ef53:	vmovaps YMMWORD PTR [rax+rdi*1],ymm7
  94ef58:	mov    rax,QWORD PTR [r14]
  94ef5b:	vmovaps ymm7,YMMWORD PTR [rsp+0xc0]
  94ef64:	mov    rdi,QWORD PTR [rsp+0x80]
  94ef6c:	vmovaps YMMWORD PTR [rax+rbx*1],ymm6
  94ef71:	mov    rax,QWORD PTR [r14]
  94ef74:	vmovaps ymm6,YMMWORD PTR [rsp+0xa0]
  94ef7d:	mov    rbx,QWORD PTR [rsp+0x70]
  94ef82:	vmovaps YMMWORD PTR [rax+r9*1],ymm7
  94ef88:	mov    rax,QWORD PTR [r14]
  94ef8b:	vmovaps ymm7,YMMWORD PTR [rsp+0xe0]
  94ef94:	vmovaps YMMWORD PTR [rax+rdi*1],ymm4
  94ef99:	mov    rdi,QWORD PTR [rsp+0x78]
  94ef9e:	mov    rax,QWORD PTR [r14]
  94efa1:	vmovaps YMMWORD PTR [rax+rdi*1],ymm6
  94efa6:	mov    rax,QWORD PTR [r14]
  94efa9:	mov    rdi,QWORD PTR [rsp+0x60]
  94efae:	vmovaps YMMWORD PTR [rax+r15*1],ymm3
  94efb4:	mov    rax,QWORD PTR [r14]
  94efb7:	vmovaps YMMWORD PTR [rax+r13*1],ymm7
  94efbd:	mov    rax,QWORD PTR [r14]
  94efc0:	vmovaps YMMWORD PTR [rax+rbx*1],ymm2
  94efc5:	mov    rax,QWORD PTR [r14]
  94efc8:	vmovaps YMMWORD PTR [rax+rdi*1],ymm1
  94efcd:	mov    rbx,QWORD PTR [rsp+0x90]
  94efd5:	add    QWORD PTR [rsp+0x98],0x8
  94efde:	mov    rax,QWORD PTR [rsp+0x98]
  94efe6:	cmp    rax,QWORD PTR [rbx+0x40]
  94efea:	jb     94ebe8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x148>
  94eff0:	vzeroupper
  94eff3:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94efff:	lea    r12,[rax-0x360]
  94f006:	mov    r13,rax
  94f009:	mov    r14,QWORD PTR [r12+0x188]
  94f011:	test   r14,r14
  94f014:	je     94f038 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x598>
  94f016:	movzx  eax,BYTE PTR [r12+0x180]
  94f01f:	lea    ebx,[rax-0x1]
  94f022:	mov    BYTE PTR [r12+0x180],bl
  94f02a:	cmp    bl,0xf
  94f02d:	jbe    94f047 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x5a7>
  94f02f:	add    WORD PTR [r14+0xfe2],0x1
  94f038:	lea    rsp,[rbp-0x28]
  94f03c:	pop    rbx
  94f03d:	pop    r12
  94f03f:	pop    r13
  94f041:	pop    r14
  94f043:	pop    r15
  94f045:	pop    rbp
  94f046:	ret
  94f047:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  94f04c:	movzx  r15d,bl
  94f050:	movsxd rcx,r15d
  94f053:	mov    rsi,rax
  94f056:	vmovq  xmm6,rax
  94f05b:	lea    rdx,[rcx+rcx*2]
  94f05f:	shl    rdx,0x3
  94f063:	sub    rsi,QWORD PTR [r12+rdx*1]
  94f067:	vpinsrq xmm0,xmm6,rsi,0x1
  94f06d:	vpsubq xmm0,xmm0,XMMWORD PTR [r13+rdx*1-0x360]
  94f077:	movzx  r13d,BYTE PTR [r12+0x181]
  94f080:	cmp    r13b,bl
  94f083:	jae    94f099 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x5f9>
  94f085:	lea    eax,[r15-0x1]
  94f089:	vmovq  rdx,xmm0
  94f08e:	cdqe
  94f090:	lea    rax,[rax+rax*2]
  94f094:	add    QWORD PTR [r12+rax*8+0x8],rdx
  94f099:	movzx  r12d,WORD PTR [r14+0xfe0]
  94f0a1:	cmp    r12w,0x7f
  94f0a6:	je     94f23c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x79c>
  94f0ac:	mov    QWORD PTR [rsp+0x120],rcx
  94f0b4:	sub    ebx,r13d
  94f0b7:	vmovdqa XMMWORD PTR [rsp+0x140],xmm0
  94f0c0:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94f0cc:	mov    rcx,QWORD PTR [rsp+0x120]
  94f0d4:	mov    DWORD PTR [rsp+0x178],0x0
  94f0df:	vmovdqa xmm0,XMMWORD PTR [rsp+0x140]
  94f0e8:	vmovdqu XMMWORD PTR [rsp+0x168],xmm0
  94f0f1:	lea    rdx,[rcx+rcx*2]
  94f0f5:	lea    rdx,[rax+rdx*8-0x360]
  94f0fd:	mov    rax,QWORD PTR [rdx]
  94f100:	mov    QWORD PTR [rsp+0x160],rax
  94f108:	movzx  eax,BYTE PTR [rdx+0x10]
  94f10c:	lea    edx,[r12+0x1]
  94f111:	vmovdqa xmm7,XMMWORD PTR [rsp+0x160]
  94f11a:	mov    WORD PTR [r14+0xfe0],dx
  94f122:	mov    ah,bl
  94f124:	mov    WORD PTR [rsp+0x17c],ax
  94f12c:	movzx  eax,r12w
  94f130:	vmovdqu xmm6,XMMWORD PTR [rsp+0x16e]
  94f139:	shl    rax,0x5
  94f13d:	add    rax,r14
  94f140:	vmovdqu XMMWORD PTR [rax],xmm7
  94f144:	vmovdqu XMMWORD PTR [rax+0xe],xmm6
  94f149:	lea    rsp,[rbp-0x28]
  94f14d:	pop    rbx
  94f14e:	pop    r12
  94f150:	pop    r13
  94f152:	pop    r14
  94f154:	pop    r15
  94f156:	pop    rbp
  94f157:	ret
  94f158:	shl    r8,0x4
  94f15c:	add    rdi,r8
  94f15f:	cmp    r9,rdi
  94f162:	je     94eb68 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc8>
  94f168:	mov    QWORD PTR [r14+0x8],rdi
  94f16c:	jmp    94eb68 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc8>
  94f171:	mov    eax,0x1
  94f176:	lock xadd DWORD PTR [rip+0x2d89e2],eax        # c27b60 <cpl::Profiling::registerRegion(char const*)::counter>
  94f17e:	mov    edx,0x1
  94f183:	mov    r13d,0x1
  94f189:	add    eax,0x2
  94f18c:	cmp    eax,0xfe
  94f191:	ja     94f1ac <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x70c>
  94f193:	lea    rdx,[rip+0x2d6e86]        # c26020 <cpl::Profiling::regions>
  94f19a:	mov    ecx,eax
  94f19c:	mov    r13d,eax
  94f19f:	lea    rdi,[rip+0xb210c]        # a012b2 <_IO_stdin_used+0x182b2>
  94f1a6:	mov    QWORD PTR [rdx+rcx*8],rdi
  94f1aa:	mov    edx,eax
  94f1ac:	xor    eax,eax
  94f1ae:	lock cmpxchg BYTE PTR [rip+0x2e3f52],dl        # c33108 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  94f1b6:	cmovne r13d,eax
  94f1ba:	jmp    94eae0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x40>
  94f1bf:	sub    rsi,rcx
  94f1c2:	mov    rdi,r14
  94f1c5:	call   939bb0 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  94f1ca:	mov    rdi,QWORD PTR [rsp+0x90]
  94f1d2:	mov    rax,QWORD PTR [rdi+0x38]
  94f1d6:	mov    rdx,QWORD PTR [rdi+0x48]
  94f1da:	jmp    94eb68 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc8>
  94f1df:	mov    QWORD PTR [rsp+0x140],rax
  94f1e7:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  94f1ec:	mov    rdx,QWORD PTR [rsp+0x140]
  94f1f4:	mov    QWORD PTR [rsp+0x168],0x0
  94f200:	mov    QWORD PTR [rsp+0x160],rax
  94f208:	lea    rax,[r12+r12*2]
  94f20c:	vmovdqa xmm7,XMMWORD PTR [rsp+0x160]
  94f215:	shl    rax,0x3
  94f219:	mov    BYTE PTR [rsp+0x170],r13b
  94f221:	vmovdqu XMMWORD PTR [rdx+rax*1-0x360],xmm7
  94f22a:	mov    BYTE PTR [rax+r15*1+0x10],r13b
  94f22f:	movzx  r12d,BYTE PTR [r15+0x180]
  94f237:	jmp    94eb0f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x6f>
  94f23c:	movzx  eax,WORD PTR [r14+0xfe2]
  94f244:	cmp    ax,0xffff
  94f248:	je     94f038 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x598>
  94f24e:	add    eax,0x1
  94f251:	mov    WORD PTR [r14+0xfe2],ax
  94f259:	jmp    94f038 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x598>
  94f25e:	endbr64
  94f262:	mov    rbx,rax
  94f265:	vzeroupper
  94f268:	call   93c790 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  94f26d:	mov    rdi,rbx
  94f270:	call   e1b10 <_Unwind_Resume@plt>

Disassembly of section .fini:
