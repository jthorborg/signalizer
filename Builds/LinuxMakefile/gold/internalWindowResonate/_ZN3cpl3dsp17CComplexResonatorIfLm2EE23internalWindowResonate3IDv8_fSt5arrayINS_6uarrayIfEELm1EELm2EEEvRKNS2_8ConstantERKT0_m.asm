; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

000000000094ea40 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long):
  94ea40:	endbr64
  94ea44:	push   rbp
  94ea45:	mov    rbp,rsp
  94ea48:	push   r15
  94ea4a:	push   r14
  94ea4c:	mov    r14,rdi
  94ea4f:	push   r13
  94ea51:	push   r12
  94ea53:	push   rbx
  94ea54:	mov    rbx,rcx
  94ea57:	and    rsp,0xffffffffffffffe0
  94ea5b:	sub    rsp,0x160
  94ea62:	mov    QWORD PTR [rsp+0x70],rsi
  94ea67:	mov    QWORD PTR [rsp+0x30],rdx
  94ea6c:	movzx  r13d,BYTE PTR [rip+0x2e3694]        # c32108 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  94ea74:	test   r13b,r13b
  94ea77:	je     94f048 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x608>
  94ea7d:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94ea89:	cmp    QWORD PTR [rax-0x1d8],0x0
  94ea91:	lea    r15,[rax-0x360]
  94ea98:	je     94eac3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x83>
  94ea9a:	movzx  r12d,BYTE PTR [r15+0x180]
  94eaa2:	cmp    r12b,0xf
  94eaa6:	jbe    94f0b3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x673>
  94eaac:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94eab8:	add    r12d,0x1
  94eabc:	mov    BYTE PTR [rax-0x1e0],r12b
  94eac3:	mov    rax,QWORD PTR [rsp+0x70]
  94eac8:	mov    r9,QWORD PTR [r14+0x8]
  94eacc:	mov    rdi,QWORD PTR [r14]
  94eacf:	mov    rdx,QWORD PTR [rax+0x48]
  94ead3:	mov    rax,QWORD PTR [rax+0x38]
  94ead7:	mov    rcx,r9
  94eada:	sub    rcx,rdi
  94eadd:	sar    rcx,0x2
  94eae1:	mov    r8,rax
  94eae4:	imul   r8,rdx
  94eae8:	lea    rsi,[r8*4+0x0]
  94eaf0:	cmp    rcx,rsi
  94eaf3:	jb     94f096 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x656>
  94eaf9:	cmp    rsi,rcx
  94eafc:	jb     94f02f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x5ef>
  94eb02:	mov    rdi,QWORD PTR [rsp+0x70]
  94eb07:	lea    rcx,[rdx+rdx*1]
  94eb0b:	imul   rax,rcx
  94eb0f:	cmp    QWORD PTR [rdi+0x40],0x0
  94eb14:	je     94eeca <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x48a>
  94eb1a:	lea    r11,[rdx*4+0x0]
  94eb22:	xor    r10d,r10d
  94eb25:	mov    QWORD PTR [rsp+0x78],0x0
  94eb2e:	lea    rsi,[r11+rax*1]
  94eb32:	shl    rax,0x2
  94eb36:	mov    QWORD PTR [rsp+0x28],r11
  94eb3b:	mov    QWORD PTR [rsp+0x10],rax
  94eb40:	lea    rdi,[rdx*8+0x0]
  94eb48:	mov    rax,rsi
  94eb4b:	shl    rdx,0x4
  94eb4f:	sub    rax,rcx
  94eb52:	mov    QWORD PTR [rsp+0x20],rdi
  94eb57:	shl    rax,0x2
  94eb5b:	mov    QWORD PTR [rsp+0x18],rdx
  94eb60:	mov    QWORD PTR [rsp+0x60],rax
  94eb65:	lea    rax,[rsi*4+0x0]
  94eb6d:	mov    QWORD PTR [rsp+0x68],rax
  94eb72:	nop    WORD PTR [rax+rax*1+0x0]
  94eb78:	mov    rdi,QWORD PTR [rsp+0x20]
  94eb7d:	mov    rcx,QWORD PTR [rsp+0x28]
  94eb82:	mov    rax,QWORD PTR [rsp+0x70]
  94eb87:	mov    rdx,QWORD PTR [rsp+0x30]
  94eb8c:	add    rdi,r10
  94eb8f:	mov    r12,QWORD PTR [rsp+0x10]
  94eb94:	lea    rsi,[rcx+rdi*1]
  94eb98:	lea    r9,[rcx+rsi*1]
  94eb9c:	mov    rcx,QWORD PTR [rsp+0x18]
  94eba1:	mov    rax,QWORD PTR [rax]
  94eba4:	lea    r13,[r10+r12*1]
  94eba8:	add    r12,r11
  94ebab:	lea    r8,[r11+rcx*1]
  94ebaf:	mov    rcx,QWORD PTR [r14]
  94ebb2:	vmovaps ymm6,YMMWORD PTR [rax+r10*1]
  94ebb8:	vmovaps ymm7,YMMWORD PTR [rax+r11*1]
  94ebbe:	vmovaps ymm9,YMMWORD PTR [rax+rdi*1]
  94ebc3:	vmovaps ymm10,YMMWORD PTR [rax+rsi*1]
  94ebc8:	lea    r15,[rcx+r10*1]
  94ebcc:	vmovaps ymm11,YMMWORD PTR [rax+r9*1]
  94ebd2:	vmovaps ymm8,YMMWORD PTR [rax+r8*1]
  94ebd8:	mov    rax,QWORD PTR [rdx]
  94ebdb:	vmovaps ymm2,YMMWORD PTR [rcx+r11*1]
  94ebe1:	vmovaps ymm3,YMMWORD PTR [rcx+rsi*1]
  94ebe6:	vmovaps YMMWORD PTR [rsp+0xa0],ymm2
  94ebef:	vmovaps ymm2,YMMWORD PTR [rcx+r9*1]
  94ebf5:	mov    QWORD PTR [rsp+0x38],rax
  94ebfa:	mov    rax,QWORD PTR [rsp+0x60]
  94ebff:	vmovaps YMMWORD PTR [rsp+0x100],ymm3
  94ec08:	vmovaps ymm3,YMMWORD PTR [rcx+r8*1]
  94ec0e:	vmovaps YMMWORD PTR [rsp+0x120],ymm2
  94ec17:	vmovaps ymm2,YMMWORD PTR [rcx+r12*1]
  94ec1d:	vmovaps YMMWORD PTR [rsp+0xc0],ymm3
  94ec26:	mov    rdx,QWORD PTR [rdx+0x10]
  94ec2a:	add    rax,r10
  94ec2d:	vmovaps ymm1,YMMWORD PTR [r15]
  94ec32:	vmovaps YMMWORD PTR [rsp+0x80],ymm2
  94ec3b:	mov    QWORD PTR [rsp+0x58],rax
  94ec40:	vmovaps ymm0,YMMWORD PTR [rcx+rdi*1]
  94ec45:	vmovaps ymm5,YMMWORD PTR [rcx+r13*1]
  94ec4b:	vmovaps ymm4,YMMWORD PTR [rcx+rax*1]
  94ec50:	mov    rax,QWORD PTR [rsp+0x60]
  94ec55:	add    rax,r11
  94ec58:	vmovaps ymm3,YMMWORD PTR [rcx+rax*1]
  94ec5d:	mov    QWORD PTR [rsp+0x50],rax
  94ec62:	mov    rax,QWORD PTR [rsp+0x68]
  94ec67:	vmovaps YMMWORD PTR [rsp+0xe0],ymm3
  94ec70:	add    rax,r10
  94ec73:	vmovaps ymm3,YMMWORD PTR [rcx+rax*1]
  94ec78:	mov    QWORD PTR [rsp+0x48],rax
  94ec7d:	mov    rax,QWORD PTR [rsp+0x68]
  94ec82:	add    rax,r11
  94ec85:	vmovaps ymm2,YMMWORD PTR [rcx+rax*1]
  94ec8a:	mov    QWORD PTR [rsp+0x40],rax
  94ec8f:	test   rbx,rbx
  94ec92:	je     94edf9 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3b9>
  94ec98:	mov    rax,QWORD PTR [rsp+0x38]
  94ec9d:	xor    r15d,r15d
  94eca0:	vmulps ymm15,ymm1,ymm7
  94eca4:	vmulps ymm1,ymm1,ymm6
  94eca8:	add    r15,0x1
  94ecac:	add    rdx,0x4
  94ecb0:	vbroadcastss ymm12,DWORD PTR [rax]
  94ecb5:	add    rax,0x4
  94ecb9:	vmovaps ymm14,YMMWORD PTR [rsp+0xa0]
  94ecc2:	vaddps ymm1,ymm1,ymm12
  94ecc7:	vmulps ymm13,ymm7,ymm14
  94eccc:	vmulps ymm14,ymm6,ymm14
  94ecd1:	vaddps ymm15,ymm15,ymm14
  94ecd6:	vmovaps ymm14,YMMWORD PTR [rsp+0x100]
  94ecdf:	vsubps ymm1,ymm1,ymm13
  94ece4:	vmulps ymm13,ymm10,ymm14
  94ece9:	vmulps ymm14,ymm9,ymm14
  94ecee:	vmovaps YMMWORD PTR [rsp+0xa0],ymm15
  94ecf7:	vmulps ymm15,ymm10,ymm0
  94ecfb:	vmulps ymm0,ymm9,ymm0
  94ecff:	vaddps ymm15,ymm15,ymm14
  94ed04:	vmovaps ymm14,YMMWORD PTR [rsp+0xc0]
  94ed0d:	vaddps ymm0,ymm0,ymm12
  94ed12:	vmovaps YMMWORD PTR [rsp+0x100],ymm15
  94ed1b:	vmulps ymm15,ymm8,YMMWORD PTR [rsp+0x120]
  94ed24:	vsubps ymm0,ymm0,ymm13
  94ed29:	vmulps ymm13,ymm8,ymm14
  94ed2e:	vmulps ymm14,ymm11,ymm14
  94ed33:	vaddps ymm14,ymm15,ymm14
  94ed38:	vmovaps ymm15,YMMWORD PTR [rsp+0x80]
  94ed41:	vmovaps YMMWORD PTR [rsp+0xc0],ymm14
  94ed4a:	vmulps ymm14,ymm11,YMMWORD PTR [rsp+0x120]
  94ed53:	vaddps ymm14,ymm14,ymm12
  94ed58:	vmulps ymm12,ymm7,ymm15
  94ed5d:	vsubps ymm13,ymm14,ymm13
  94ed62:	vmulps ymm14,ymm6,ymm5
  94ed66:	vmulps ymm5,ymm7,ymm5
  94ed6a:	vsubps ymm12,ymm14,ymm12
  94ed6f:	vmulps ymm14,ymm6,ymm15
  94ed74:	vmovaps ymm15,YMMWORD PTR [rsp+0xe0]
  94ed7d:	vaddps ymm5,ymm5,ymm14
  94ed82:	vmovaps YMMWORD PTR [rsp+0x120],ymm13
  94ed8b:	vbroadcastss ymm13,DWORD PTR [rdx-0x4]
  94ed91:	vmulps ymm14,ymm10,ymm15
  94ed96:	vmovaps YMMWORD PTR [rsp+0x80],ymm5
  94ed9f:	vaddps ymm5,ymm13,ymm12
  94eda4:	vmulps ymm12,ymm4,ymm9
  94eda9:	vmulps ymm4,ymm4,ymm10
  94edae:	vsubps ymm12,ymm12,ymm14
  94edb3:	vmulps ymm14,ymm9,ymm15
  94edb8:	vaddps ymm4,ymm4,ymm14
  94edbd:	vmulps ymm14,ymm3,ymm11
  94edc2:	vmulps ymm3,ymm8,ymm3
  94edc6:	vmovaps YMMWORD PTR [rsp+0xe0],ymm4
  94edcf:	vaddps ymm4,ymm13,ymm12
  94edd4:	vmulps ymm12,ymm2,ymm8
  94edd9:	vmulps ymm2,ymm2,ymm11
  94edde:	vsubps ymm12,ymm14,ymm12
  94ede3:	vaddps ymm2,ymm3,ymm2
  94ede7:	vaddps ymm3,ymm13,ymm12
  94edec:	cmp    rbx,r15
  94edef:	jne    94eca0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x260>
  94edf5:	lea    r15,[rcx+r10*1]
  94edf9:	vmovaps YMMWORD PTR [r15],ymm1
  94edfe:	mov    rax,QWORD PTR [r14]
  94ee01:	add    r10,0x20
  94ee05:	vmovaps ymm6,YMMWORD PTR [rsp+0xa0]
  94ee0e:	vmovaps ymm7,YMMWORD PTR [rsp+0x100]
  94ee17:	vmovaps YMMWORD PTR [rax+r11*1],ymm6
  94ee1d:	mov    rax,QWORD PTR [r14]
  94ee20:	add    r11,0x20
  94ee24:	vmovaps ymm6,YMMWORD PTR [rsp+0x120]
  94ee2d:	vmovaps YMMWORD PTR [rax+rdi*1],ymm0
  94ee32:	mov    rax,QWORD PTR [r14]
  94ee35:	mov    rdi,QWORD PTR [rsp+0x50]
  94ee3a:	vmovaps YMMWORD PTR [rax+rsi*1],ymm7
  94ee3f:	mov    rax,QWORD PTR [r14]
  94ee42:	vmovaps ymm7,YMMWORD PTR [rsp+0xc0]
  94ee4b:	mov    rsi,QWORD PTR [rsp+0x58]
  94ee50:	vmovaps YMMWORD PTR [rax+r9*1],ymm6
  94ee56:	mov    rax,QWORD PTR [r14]
  94ee59:	vmovaps ymm6,YMMWORD PTR [rsp+0x80]
  94ee62:	vmovaps YMMWORD PTR [rax+r8*1],ymm7
  94ee68:	mov    rax,QWORD PTR [r14]
  94ee6b:	vmovaps ymm7,YMMWORD PTR [rsp+0xe0]
  94ee74:	vmovaps YMMWORD PTR [rax+r13*1],ymm5
  94ee7a:	mov    rax,QWORD PTR [r14]
  94ee7d:	vmovaps YMMWORD PTR [rax+r12*1],ymm6
  94ee83:	mov    rax,QWORD PTR [r14]
  94ee86:	vmovaps YMMWORD PTR [rax+rsi*1],ymm4
  94ee8b:	mov    rax,QWORD PTR [r14]
  94ee8e:	mov    rsi,QWORD PTR [rsp+0x48]
  94ee93:	vmovaps YMMWORD PTR [rax+rdi*1],ymm7
  94ee98:	mov    rax,QWORD PTR [r14]
  94ee9b:	mov    rdi,QWORD PTR [rsp+0x40]
  94eea0:	vmovaps YMMWORD PTR [rax+rsi*1],ymm3
  94eea5:	mov    rax,QWORD PTR [r14]
  94eea8:	vmovaps YMMWORD PTR [rax+rdi*1],ymm2
  94eead:	mov    rsi,QWORD PTR [rsp+0x70]
  94eeb2:	add    QWORD PTR [rsp+0x78],0x8
  94eeb8:	mov    rax,QWORD PTR [rsp+0x78]
  94eebd:	cmp    rax,QWORD PTR [rsi+0x40]
  94eec1:	jb     94eb78 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x138>
  94eec7:	vzeroupper
  94eeca:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94eed6:	lea    r12,[rax-0x360]
  94eedd:	mov    r13,rax
  94eee0:	mov    r14,QWORD PTR [r12+0x188]
  94eee8:	test   r14,r14
  94eeeb:	je     94ef0f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4cf>
  94eeed:	movzx  eax,BYTE PTR [r12+0x180]
  94eef6:	lea    ebx,[rax-0x1]
  94eef9:	mov    BYTE PTR [r12+0x180],bl
  94ef01:	cmp    bl,0xf
  94ef04:	jbe    94ef1e <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4de>
  94ef06:	add    WORD PTR [r14+0xfe2],0x1
  94ef0f:	lea    rsp,[rbp-0x28]
  94ef13:	pop    rbx
  94ef14:	pop    r12
  94ef16:	pop    r13
  94ef18:	pop    r14
  94ef1a:	pop    r15
  94ef1c:	pop    rbp
  94ef1d:	ret
  94ef1e:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  94ef23:	movzx  r15d,bl
  94ef27:	movsxd rcx,r15d
  94ef2a:	mov    rsi,rax
  94ef2d:	vmovq  xmm6,rax
  94ef32:	lea    rdx,[rcx+rcx*2]
  94ef36:	shl    rdx,0x3
  94ef3a:	sub    rsi,QWORD PTR [r12+rdx*1]
  94ef3e:	vpinsrq xmm0,xmm6,rsi,0x1
  94ef44:	vpsubq xmm0,xmm0,XMMWORD PTR [r13+rdx*1-0x360]
  94ef4e:	movzx  r13d,BYTE PTR [r12+0x181]
  94ef57:	cmp    r13b,bl
  94ef5a:	jae    94ef70 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x530>
  94ef5c:	lea    eax,[r15-0x1]
  94ef60:	vmovq  rdx,xmm0
  94ef65:	cdqe
  94ef67:	lea    rax,[rax+rax*2]
  94ef6b:	add    QWORD PTR [r12+rax*8+0x8],rdx
  94ef70:	movzx  r12d,WORD PTR [r14+0xfe0]
  94ef78:	cmp    r12w,0x7f
  94ef7d:	je     94f110 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x6d0>
  94ef83:	mov    QWORD PTR [rsp+0x100],rcx
  94ef8b:	sub    ebx,r13d
  94ef8e:	vmovdqa XMMWORD PTR [rsp+0x120],xmm0
  94ef97:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94efa3:	mov    rcx,QWORD PTR [rsp+0x100]
  94efab:	mov    DWORD PTR [rsp+0x158],0x0
  94efb6:	vmovdqa xmm0,XMMWORD PTR [rsp+0x120]
  94efbf:	vmovdqu XMMWORD PTR [rsp+0x148],xmm0
  94efc8:	lea    rdx,[rcx+rcx*2]
  94efcc:	lea    rdx,[rax+rdx*8-0x360]
  94efd4:	mov    rax,QWORD PTR [rdx]
  94efd7:	mov    QWORD PTR [rsp+0x140],rax
  94efdf:	movzx  eax,BYTE PTR [rdx+0x10]
  94efe3:	lea    edx,[r12+0x1]
  94efe8:	vmovdqa xmm7,XMMWORD PTR [rsp+0x140]
  94eff1:	mov    WORD PTR [r14+0xfe0],dx
  94eff9:	mov    ah,bl
  94effb:	mov    WORD PTR [rsp+0x15c],ax
  94f003:	movzx  eax,r12w
  94f007:	vmovdqu xmm6,XMMWORD PTR [rsp+0x14e]
  94f010:	shl    rax,0x5
  94f014:	add    rax,r14
  94f017:	vmovdqu XMMWORD PTR [rax],xmm7
  94f01b:	vmovdqu XMMWORD PTR [rax+0xe],xmm6
  94f020:	lea    rsp,[rbp-0x28]
  94f024:	pop    rbx
  94f025:	pop    r12
  94f027:	pop    r13
  94f029:	pop    r14
  94f02b:	pop    r15
  94f02d:	pop    rbp
  94f02e:	ret
  94f02f:	shl    r8,0x4
  94f033:	add    rdi,r8
  94f036:	cmp    r9,rdi
  94f039:	je     94eb02 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc2>
  94f03f:	mov    QWORD PTR [r14+0x8],rdi
  94f043:	jmp    94eb02 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc2>
  94f048:	mov    eax,0x1
  94f04d:	lock xadd DWORD PTR [rip+0x2d7b0b],eax        # c26b60 <cpl::Profiling::registerRegion(char const*)::counter>
  94f055:	mov    edx,0x1
  94f05a:	mov    r13d,0x1
  94f060:	add    eax,0x2
  94f063:	cmp    eax,0xfe
  94f068:	ja     94f083 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x643>
  94f06a:	lea    rdx,[rip+0x2d5faf]        # c25020 <cpl::Profiling::regions>
  94f071:	mov    ecx,eax
  94f073:	mov    r13d,eax
  94f076:	lea    rsi,[rip+0xb1235]        # a002b2 <_IO_stdin_used+0x182b2>
  94f07d:	mov    QWORD PTR [rdx+rcx*8],rsi
  94f081:	mov    edx,eax
  94f083:	xor    eax,eax
  94f085:	lock cmpxchg BYTE PTR [rip+0x2e307b],dl        # c32108 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  94f08d:	cmovne r13d,eax
  94f091:	jmp    94ea7d <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3d>
  94f096:	sub    rsi,rcx
  94f099:	mov    rdi,r14
  94f09c:	call   939b50 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  94f0a1:	mov    rdi,QWORD PTR [rsp+0x70]
  94f0a6:	mov    rax,QWORD PTR [rdi+0x38]
  94f0aa:	mov    rdx,QWORD PTR [rdi+0x48]
  94f0ae:	jmp    94eb02 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc2>
  94f0b3:	mov    QWORD PTR [rsp+0x120],rax
  94f0bb:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  94f0c0:	mov    rdx,QWORD PTR [rsp+0x120]
  94f0c8:	mov    QWORD PTR [rsp+0x148],0x0
  94f0d4:	mov    QWORD PTR [rsp+0x140],rax
  94f0dc:	lea    rax,[r12+r12*2]
  94f0e0:	vmovdqa xmm7,XMMWORD PTR [rsp+0x140]
  94f0e9:	shl    rax,0x3
  94f0ed:	mov    BYTE PTR [rsp+0x150],r13b
  94f0f5:	vmovdqu XMMWORD PTR [rdx+rax*1-0x360],xmm7
  94f0fe:	mov    BYTE PTR [rax+r15*1+0x10],r13b
  94f103:	movzx  r12d,BYTE PTR [r15+0x180]
  94f10b:	jmp    94eaac <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x6c>
  94f110:	movzx  eax,WORD PTR [r14+0xfe2]
  94f118:	cmp    ax,0xffff
  94f11c:	je     94ef0f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4cf>
  94f122:	add    eax,0x1
  94f125:	mov    WORD PTR [r14+0xfe2],ax
  94f12d:	jmp    94ef0f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4cf>
  94f132:	endbr64
  94f136:	mov    rbx,rax
  94f139:	vzeroupper
  94f13c:	call   93c730 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  94f141:	mov    rdi,rbx
  94f144:	call   e1b10 <_Unwind_Resume@plt>

Disassembly of section .fini:
