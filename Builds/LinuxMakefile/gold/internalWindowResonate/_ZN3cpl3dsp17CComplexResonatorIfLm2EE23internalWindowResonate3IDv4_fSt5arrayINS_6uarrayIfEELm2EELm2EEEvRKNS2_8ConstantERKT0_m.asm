; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000944990 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long):
  944990:	endbr64
  944994:	push   r15
  944996:	push   r14
  944998:	push   r13
  94499a:	push   r12
  94499c:	mov    r12,rdi
  94499f:	push   rbp
  9449a0:	push   rbx
  9449a1:	sub    rsp,0x298
  9449a8:	mov    QWORD PTR [rsp+0x80],rsi
  9449b0:	mov    QWORD PTR [rsp+0x88],rdx
  9449b8:	mov    QWORD PTR [rsp+0xc8],rcx
  9449c0:	mov    rax,QWORD PTR fs:0x28
  9449c9:	mov    QWORD PTR [rsp+0x288],rax
  9449d1:	xor    eax,eax
  9449d3:	movzx  r13d,BYTE PTR [rip+0x2da70d]        # c1f0e8 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  9449db:	test   r13b,r13b
  9449de:	je     945237 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x8a7>
  9449e4:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  9449f0:	cmp    QWORD PTR [rax-0x1d8],0x0
  9449f8:	mov    rbx,rax
  9449fb:	lea    r14,[rax-0x360]
  944a02:	je     944a2c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x9c>
  944a04:	movzx  ebp,BYTE PTR [r14+0x180]
  944a0c:	cmp    bpl,0xf
  944a10:	jbe    9452a5 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x915>
  944a16:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  944a22:	add    ebp,0x1
  944a25:	mov    BYTE PTR [rax-0x1e0],bpl
  944a2c:	mov    rsi,QWORD PTR [rsp+0x80]
  944a34:	mov    r9,QWORD PTR [r12+0x8]
  944a39:	mov    rdi,QWORD PTR [r12]
  944a3d:	mov    rdx,QWORD PTR [rsi+0x38]
  944a41:	mov    rax,QWORD PTR [rsi+0x48]
  944a45:	mov    rcx,r9
  944a48:	sub    rcx,rdi
  944a4b:	mov    r8,rdx
  944a4e:	sar    rcx,0x2
  944a52:	imul   r8,rax
  944a56:	lea    rsi,[r8*4+0x0]
  944a5e:	cmp    rcx,rsi
  944a61:	jb     945285 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x8f5>
  944a67:	cmp    rsi,rcx
  944a6a:	jb     94521d <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x88d>
  944a70:	mov    rcx,QWORD PTR [rsp+0x80]
  944a78:	lea    rsi,[rax+rax*1]
  944a7c:	imul   rdx,rsi
  944a80:	cmp    QWORD PTR [rcx+0x40],0x0
  944a85:	je     9450bf <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x72f>
  944a8b:	lea    r8,[rsi+rax*1]
  944a8f:	lea    rcx,[rax*4+0x0]
  944a97:	xor    r14d,r14d
  944a9a:	mov    QWORD PTR [rsp+0x70],0x0
  944aa3:	lea    rbx,[r8*4+0x0]
  944aab:	lea    rdi,[rcx+rdx*1]
  944aaf:	mov    QWORD PTR [rsp+0x198],rcx
  944ab7:	mov    r13,rcx
  944aba:	mov    QWORD PTR [rsp+0x60],rbx
  944abf:	mov    rbx,rax
  944ac2:	add    rax,rcx
  944ac5:	shl    rax,0x2
  944ac9:	shl    rbx,0x4
  944acd:	mov    QWORD PTR [rsp+0x50],rax
  944ad2:	lea    rax,[rdx*4+0x0]
  944ada:	mov    QWORD PTR [rsp+0x68],rax
  944adf:	mov    rax,rdi
  944ae2:	sub    rax,rsi
  944ae5:	mov    QWORD PTR [rsp+0x58],rbx
  944aea:	shl    rax,0x2
  944aee:	mov    QWORD PTR [rsp+0x98],rax
  944af6:	lea    rax,[rdi*4+0x0]
  944afe:	mov    QWORD PTR [rsp+0x90],rax
  944b06:	lea    rax,[rsp+0x1c0]
  944b0e:	mov    QWORD PTR [rsp+0xc0],rax
  944b16:	lea    rax,[rsp+0x1e0]
  944b1e:	mov    QWORD PTR [rsp+0xa0],rax
  944b26:	lea    rax,[rsp+0x200]
  944b2e:	mov    QWORD PTR [rsp+0xa8],rax
  944b36:	lea    rax,[rsp+0x220]
  944b3e:	mov    QWORD PTR [rsp+0xb0],rax
  944b46:	lea    rax,[rsp+0x240]
  944b4e:	mov    QWORD PTR [rsp+0xb8],rax
  944b56:	lea    rax,[rsp+0x260]
  944b5e:	mov    QWORD PTR [rsp+0x190],rax
  944b66:	cs nop WORD PTR [rax+rax*1+0x0]
  944b70:	mov    rax,QWORD PTR [rsp+0x80]
  944b78:	mov    rbx,QWORD PTR [rax]
  944b7b:	lea    rdi,[rbx+r14*1]
  944b7f:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944b84:	lea    rdi,[rbx+r13*1]
  944b88:	movaps XMMWORD PTR [rsp],xmm0
  944b8c:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944b91:	mov    rax,QWORD PTR [rsp+0x198]
  944b99:	movaps XMMWORD PTR [rsp+0x10],xmm0
  944b9e:	lea    r15,[rax+r13*1]
  944ba2:	lea    rdi,[rbx+r15*1]
  944ba6:	mov    QWORD PTR [rsp+0x78],r15
  944bab:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944bb0:	mov    rax,QWORD PTR [rsp+0x60]
  944bb5:	movaps XMMWORD PTR [rsp+0x20],xmm0
  944bba:	lea    rdi,[rbx+rax*1]
  944bbe:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944bc3:	mov    rsi,QWORD PTR [rsp+0x58]
  944bc8:	movaps XMMWORD PTR [rsp+0x30],xmm0
  944bcd:	lea    rdi,[rbx+rsi*1]
  944bd1:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944bd6:	mov    rcx,QWORD PTR [rsp+0x50]
  944bdb:	movaps XMMWORD PTR [rsp+0x40],xmm0
  944be0:	lea    rdi,[rbx+rcx*1]
  944be4:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944be9:	mov    rbx,QWORD PTR [r12]
  944bed:	mov    rdx,QWORD PTR [rsp+0x88]
  944bf5:	movaps XMMWORD PTR [rsp+0x180],xmm0
  944bfd:	lea    rdi,[rbx+r14*1]
  944c01:	mov    rbp,QWORD PTR [rdx]
  944c04:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944c09:	lea    rdi,[rbx+r13*1]
  944c0d:	movaps XMMWORD PTR [rsp+0x1c0],xmm0
  944c15:	movaps XMMWORD PTR [rsp+0x170],xmm0
  944c1d:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944c22:	lea    rdi,[rbx+r15*1]
  944c26:	movaps XMMWORD PTR [rsp+0x1e0],xmm0
  944c2e:	movaps XMMWORD PTR [rsp+0x160],xmm0
  944c36:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944c3b:	mov    rax,QWORD PTR [rsp+0x60]
  944c40:	movaps XMMWORD PTR [rsp+0x200],xmm0
  944c48:	lea    rdi,[rbx+rax*1]
  944c4c:	movaps XMMWORD PTR [rsp+0x150],xmm0
  944c54:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944c59:	mov    rsi,QWORD PTR [rsp+0x58]
  944c5e:	movaps XMMWORD PTR [rsp+0x220],xmm0
  944c66:	lea    rdi,[rbx+rsi*1]
  944c6a:	movaps XMMWORD PTR [rsp+0x140],xmm0
  944c72:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944c77:	mov    rcx,QWORD PTR [rsp+0x50]
  944c7c:	movaps XMMWORD PTR [rsp+0x240],xmm0
  944c84:	lea    rdi,[rbx+rcx*1]
  944c88:	movaps XMMWORD PTR [rsp+0x130],xmm0
  944c90:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944c95:	mov    rax,QWORD PTR [rsp+0x68]
  944c9a:	mov    rdx,QWORD PTR [rsp+0x88]
  944ca2:	movaps XMMWORD PTR [rsp+0x260],xmm0
  944caa:	lea    rdi,[rax+r14*1]
  944cae:	mov    r15,QWORD PTR [rdx+0x10]
  944cb2:	movaps XMMWORD PTR [rsp+0x120],xmm0
  944cba:	add    rdi,rbx
  944cbd:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944cc2:	mov    rax,QWORD PTR [rsp+0x68]
  944cc7:	movaps XMMWORD PTR [rsp+0x1d0],xmm0
  944ccf:	lea    rdi,[rax+r13*1]
  944cd3:	movaps XMMWORD PTR [rsp+0x110],xmm0
  944cdb:	add    rdi,rbx
  944cde:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944ce3:	movaps XMMWORD PTR [rsp+0x1f0],xmm0
  944ceb:	movaps XMMWORD PTR [rsp+0x100],xmm0
  944cf3:	mov    rax,QWORD PTR [rsp+0x98]
  944cfb:	lea    rdi,[rax+r14*1]
  944cff:	add    rdi,rbx
  944d02:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944d07:	mov    rax,QWORD PTR [rsp+0x98]
  944d0f:	movaps XMMWORD PTR [rsp+0x210],xmm0
  944d17:	lea    rdi,[rax+r13*1]
  944d1b:	movaps XMMWORD PTR [rsp+0xf0],xmm0
  944d23:	add    rdi,rbx
  944d26:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944d2b:	mov    rax,QWORD PTR [rsp+0x90]
  944d33:	movaps XMMWORD PTR [rsp+0x230],xmm0
  944d3b:	lea    rdi,[rax+r14*1]
  944d3f:	movaps XMMWORD PTR [rsp+0xe0],xmm0
  944d47:	add    rdi,rbx
  944d4a:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944d4f:	mov    rax,QWORD PTR [rsp+0x90]
  944d57:	movaps XMMWORD PTR [rsp+0x250],xmm0
  944d5f:	lea    rdi,[rax+r13*1]
  944d63:	movaps XMMWORD PTR [rsp+0xd0],xmm0
  944d6b:	add    rdi,rbx
  944d6e:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944d73:	xor    eax,eax
  944d75:	cmp    QWORD PTR [rsp+0xc8],0x0
  944d7e:	movaps xmm1,XMMWORD PTR [rsp+0xd0]
  944d86:	movaps xmm4,XMMWORD PTR [rsp+0xe0]
  944d8e:	movaps xmm3,XMMWORD PTR [rsp+0xf0]
  944d96:	movaps xmm2,xmm0
  944d99:	movaps XMMWORD PTR [rsp+0x270],xmm0
  944da1:	movaps xmm6,XMMWORD PTR [rsp+0x100]
  944da9:	movaps xmm5,XMMWORD PTR [rsp+0x110]
  944db1:	movaps xmm8,XMMWORD PTR [rsp+0x120]
  944dba:	movaps xmm7,XMMWORD PTR [rsp+0x130]
  944dc2:	movaps xmm10,XMMWORD PTR [rsp+0x140]
  944dcb:	movaps xmm9,XMMWORD PTR [rsp+0x150]
  944dd4:	movaps xmm12,XMMWORD PTR [rsp+0x160]
  944ddd:	movaps xmm11,XMMWORD PTR [rsp+0x170]
  944de6:	je     944f90 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x600>
  944dec:	movaps xmm13,XMMWORD PTR [rsp+0x180]
  944df5:	mov    rdx,QWORD PTR [rsp+0xc8]
  944dfd:	nop    DWORD PTR [rax]
  944e00:	movss  xmm0,DWORD PTR [rbp+rax*4+0x0]
  944e06:	call   9205e0 <float __vector(4) cpl::simd::broadcast<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*) [clone .isra.0]>
  944e0b:	movaps xmm14,XMMWORD PTR [rsp]
  944e10:	movaps xmm15,XMMWORD PTR [rsp+0x10]
  944e16:	mulps  xmm14,xmm11
  944e1a:	mulps  xmm15,xmm12
  944e1e:	mulps  xmm11,XMMWORD PTR [rsp+0x10]
  944e24:	mulps  xmm12,XMMWORD PTR [rsp]
  944e29:	subps  xmm14,xmm15
  944e2d:	movaps xmm15,XMMWORD PTR [rsp+0x30]
  944e33:	addps  xmm12,xmm11
  944e37:	mulps  xmm15,xmm10
  944e3b:	mulps  xmm10,XMMWORD PTR [rsp+0x20]
  944e41:	movaps xmm11,xmm14
  944e45:	movaps xmm14,XMMWORD PTR [rsp+0x20]
  944e4b:	addps  xmm11,xmm0
  944e4f:	mulps  xmm14,xmm9
  944e53:	mulps  xmm9,XMMWORD PTR [rsp+0x30]
  944e59:	subps  xmm14,xmm15
  944e5d:	movaps xmm15,xmm8
  944e61:	addps  xmm10,xmm9
  944e65:	mulps  xmm15,xmm13
  944e69:	mulps  xmm8,XMMWORD PTR [rsp+0x40]
  944e6f:	movaps xmm9,xmm14
  944e73:	movaps xmm14,XMMWORD PTR [rsp+0x40]
  944e79:	addps  xmm9,xmm0
  944e7d:	mulps  xmm14,xmm7
  944e81:	mulps  xmm7,xmm13
  944e85:	subps  xmm14,xmm15
  944e89:	addps  xmm8,xmm7
  944e8d:	movaps xmm7,xmm0
  944e90:	movss  xmm0,DWORD PTR [r15+rax*4]
  944e96:	add    rax,0x1
  944e9a:	call   9205e0 <float __vector(4) cpl::simd::broadcast<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*) [clone .isra.0]>
  944e9f:	addps  xmm7,xmm14
  944ea3:	movaps xmm15,XMMWORD PTR [rsp+0x10]
  944ea9:	movaps xmm14,XMMWORD PTR [rsp]
  944eae:	mulps  xmm15,xmm6
  944eb2:	mulps  xmm14,xmm5
  944eb6:	mulps  xmm6,XMMWORD PTR [rsp]
  944eba:	mulps  xmm5,XMMWORD PTR [rsp+0x10]
  944ebf:	subps  xmm14,xmm15
  944ec3:	movaps xmm15,XMMWORD PTR [rsp+0x30]
  944ec9:	addps  xmm6,xmm5
  944ecc:	mulps  xmm15,xmm4
  944ed0:	mulps  xmm4,XMMWORD PTR [rsp+0x20]
  944ed5:	movaps xmm5,xmm14
  944ed9:	movaps xmm14,XMMWORD PTR [rsp+0x20]
  944edf:	addps  xmm5,xmm0
  944ee2:	mulps  xmm14,xmm3
  944ee6:	mulps  xmm3,XMMWORD PTR [rsp+0x30]
  944eeb:	subps  xmm14,xmm15
  944eef:	movaps xmm15,xmm13
  944ef3:	addps  xmm4,xmm3
  944ef6:	mulps  xmm15,xmm2
  944efa:	mulps  xmm2,XMMWORD PTR [rsp+0x40]
  944eff:	movaps xmm3,xmm14
  944f03:	movaps xmm14,XMMWORD PTR [rsp+0x40]
  944f09:	addps  xmm3,xmm0
  944f0c:	mulps  xmm14,xmm1
  944f10:	mulps  xmm1,xmm13
  944f14:	subps  xmm14,xmm15
  944f18:	addps  xmm2,xmm1
  944f1b:	movaps xmm1,xmm0
  944f1e:	addps  xmm1,xmm14
  944f22:	cmp    rdx,rax
  944f25:	jne    944e00 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x470>
  944f2b:	movaps XMMWORD PTR [rsp+0x1e0],xmm12
  944f34:	movaps XMMWORD PTR [rsp+0x1c0],xmm11
  944f3d:	movaps XMMWORD PTR [rsp+0x220],xmm10
  944f46:	movaps XMMWORD PTR [rsp+0x200],xmm9
  944f4f:	movaps XMMWORD PTR [rsp+0x260],xmm8
  944f58:	movaps XMMWORD PTR [rsp+0x240],xmm7
  944f60:	movaps XMMWORD PTR [rsp+0x1f0],xmm6
  944f68:	movaps XMMWORD PTR [rsp+0x1d0],xmm5
  944f70:	movaps XMMWORD PTR [rsp+0x230],xmm4
  944f78:	movaps XMMWORD PTR [rsp+0x210],xmm3
  944f80:	movaps XMMWORD PTR [rsp+0x270],xmm2
  944f88:	movaps XMMWORD PTR [rsp+0x250],xmm1
  944f90:	xor    r15d,r15d
  944f93:	mov    rdi,rbx
  944f96:	xor    ebp,ebp
  944f98:	mov    rbx,r15
  944f9b:	mov    r15,QWORD PTR [rsp+0x190]
  944fa3:	mov    rax,QWORD PTR [rsp+0xc0]
  944fab:	lea    rdx,[r14+rbp*1]
  944faf:	add    rdi,rdx
  944fb2:	movaps xmm2,XMMWORD PTR [rax+rbx*1]
  944fb6:	movaps xmm0,xmm2
  944fb9:	movaps XMMWORD PTR [rsp],xmm2
  944fbd:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  944fc2:	mov    rax,QWORD PTR [rsp+0xa0]
  944fca:	lea    rdi,[r13+rbp*1+0x0]
  944fcf:	add    rdi,QWORD PTR [r12]
  944fd3:	movaps xmm4,XMMWORD PTR [rax+rbx*1]
  944fd7:	movaps xmm0,xmm4
  944fda:	movaps XMMWORD PTR [rsp],xmm4
  944fde:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  944fe3:	mov    rax,QWORD PTR [rsp+0xa8]
  944feb:	movaps xmm6,XMMWORD PTR [rax+rbx*1]
  944fef:	mov    rax,QWORD PTR [rsp+0x78]
  944ff4:	movaps xmm0,xmm6
  944ff7:	lea    rdi,[rax+rbp*1]
  944ffb:	add    rdi,QWORD PTR [r12]
  944fff:	movaps XMMWORD PTR [rsp],xmm6
  945003:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  945008:	mov    rax,QWORD PTR [rsp+0xb0]
  945010:	movaps xmm2,XMMWORD PTR [rax+rbx*1]
  945014:	mov    rax,QWORD PTR [rsp+0x60]
  945019:	movaps xmm0,xmm2
  94501c:	lea    rdi,[rax+rbp*1]
  945020:	add    rdi,QWORD PTR [r12]
  945024:	movaps XMMWORD PTR [rsp],xmm2
  945028:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  94502d:	mov    rax,QWORD PTR [rsp+0xb8]
  945035:	movaps xmm4,XMMWORD PTR [rax+rbx*1]
  945039:	mov    rax,QWORD PTR [rsp+0x58]
  94503e:	movaps xmm0,xmm4
  945041:	lea    rdi,[rax+rbp*1]
  945045:	add    rdi,QWORD PTR [r12]
  945049:	movaps XMMWORD PTR [rsp],xmm4
  94504d:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  945052:	mov    rax,QWORD PTR [rsp+0x50]
  945057:	movaps xmm0,XMMWORD PTR [r15+rbx*1]
  94505c:	add    rbx,0x10
  945060:	lea    rdi,[rax+rbp*1]
  945064:	add    rdi,QWORD PTR [r12]
  945068:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  94506d:	cmp    rbx,0x20
  945071:	je     945088 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x6f8>
  945073:	mov    rax,QWORD PTR [rsp+0x68]
  945078:	mov    rdi,QWORD PTR [r12]
  94507c:	add    rbp,rax
  94507f:	jmp    944fa3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x613>
  945084:	nop    DWORD PTR [rax+0x0]
  945088:	add    QWORD PTR [rsp+0x70],0x4
  94508e:	add    r14,0x10
  945092:	mov    rax,QWORD PTR [rsp+0x70]
  945097:	add    r13,0x10
  94509b:	mov    rcx,QWORD PTR [rsp+0x80]
  9450a3:	add    QWORD PTR [rsp+0x60],0x10
  9450a9:	add    QWORD PTR [rsp+0x58],0x10
  9450af:	add    QWORD PTR [rsp+0x50],0x10
  9450b5:	cmp    rax,QWORD PTR [rcx+0x40]
  9450b9:	jb     944b70 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x1e0>
  9450bf:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  9450cb:	lea    rbp,[rax-0x360]
  9450d2:	mov    r12,rax
  9450d5:	mov    r13,QWORD PTR [rbp+0x188]
  9450dc:	test   r13,r13
  9450df:	je     9450ff <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x76f>
  9450e1:	movzx  eax,BYTE PTR [rbp+0x180]
  9450e8:	lea    ebx,[rax-0x1]
  9450eb:	mov    BYTE PTR [rbp+0x180],bl
  9450f1:	cmp    bl,0xf
  9450f4:	jbe    945128 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x798>
  9450f6:	add    WORD PTR [r13+0xfe2],0x1
  9450ff:	mov    rax,QWORD PTR [rsp+0x288]
  945107:	sub    rax,QWORD PTR fs:0x28
  945110:	jne    94532c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x99c>
  945116:	add    rsp,0x298
  94511d:	pop    rbx
  94511e:	pop    rbp
  94511f:	pop    r12
  945121:	pop    r13
  945123:	pop    r14
  945125:	pop    r15
  945127:	ret
  945128:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  94512d:	movzx  r15d,bl
  945131:	lea    rdx,[r15+r15*2]
  945135:	mov    rcx,rax
  945138:	movq   xmm0,rax
  94513d:	shl    rdx,0x3
  945141:	movdqu xmm1,XMMWORD PTR [r12+rdx*1-0x360]
  94514b:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  945150:	movq   xmm6,rcx
  945155:	movzx  r12d,BYTE PTR [rbp+0x181]
  94515d:	punpcklqdq xmm0,xmm6
  945161:	psubq  xmm0,xmm1
  945165:	cmp    r12b,bl
  945168:	jae    945188 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x7f8>
  94516a:	lea    eax,[r15-0x1]
  94516e:	movdqa xmm1,xmm0
  945172:	cdqe
  945174:	lea    rax,[rax+rax*2]
  945178:	movq   xmm2,QWORD PTR [rbp+rax*8+0x8]
  94517e:	paddq  xmm1,xmm2
  945182:	movq   QWORD PTR [rbp+rax*8+0x8],xmm1
  945188:	movzx  ebp,WORD PTR [r13+0xfe0]
  945190:	cmp    bp,0x7f
  945194:	je     9452f2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x962>
  94519a:	movaps XMMWORD PTR [rsp],xmm0
  94519e:	sub    ebx,r12d
  9451a1:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  9451ad:	lea    rdx,[r15+r15*2]
  9451b1:	movdqa xmm0,XMMWORD PTR [rsp]
  9451b6:	mov    DWORD PTR [rsp+0x1b8],0x0
  9451c1:	movups XMMWORD PTR [rsp+0x1a8],xmm0
  9451c9:	lea    rdx,[rax+rdx*8-0x360]
  9451d1:	mov    rax,QWORD PTR [rdx]
  9451d4:	mov    QWORD PTR [rsp+0x1a0],rax
  9451dc:	movzx  eax,BYTE PTR [rdx+0x10]
  9451e0:	lea    edx,[rbp+0x1]
  9451e3:	movdqa xmm6,XMMWORD PTR [rsp+0x1a0]
  9451ec:	mov    WORD PTR [r13+0xfe0],dx
  9451f4:	mov    ah,bl
  9451f6:	mov    WORD PTR [rsp+0x1bc],ax
  9451fe:	movzx  eax,bp
  945201:	shl    rax,0x5
  945205:	add    rax,r13
  945208:	movups XMMWORD PTR [rax],xmm6
  94520b:	movdqu xmm6,XMMWORD PTR [rsp+0x1ae]
  945214:	movups XMMWORD PTR [rax+0xe],xmm6
  945218:	jmp    9450ff <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x76f>
  94521d:	shl    r8,0x4
  945221:	add    rdi,r8
  945224:	cmp    r9,rdi
  945227:	je     944a70 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xe0>
  94522d:	mov    QWORD PTR [r12+0x8],rdi
  945232:	jmp    944a70 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xe0>
  945237:	mov    eax,0x1
  94523c:	lock xadd DWORD PTR [rip+0x2ce91c],eax        # c13b60 <cpl::Profiling::registerRegion(char const*)::counter>
  945244:	mov    edx,0x1
  945249:	mov    r13d,0x1
  94524f:	add    eax,0x2
  945252:	cmp    eax,0xfe
  945257:	ja     945272 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x8e2>
  945259:	lea    rdx,[rip+0x2ccdc0]        # c12020 <cpl::Profiling::regions>
  945260:	mov    ecx,eax
  945262:	lea    rsi,[rip+0xa2edd]        # 9e8146 <_IO_stdin_used+0x18146>
  945269:	mov    r13d,eax
  94526c:	mov    QWORD PTR [rdx+rcx*8],rsi
  945270:	mov    edx,eax
  945272:	xor    eax,eax
  945274:	lock cmpxchg BYTE PTR [rip+0x2d9e6c],dl        # c1f0e8 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  94527c:	cmovne r13d,eax
  945280:	jmp    9449e4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x54>
  945285:	sub    rsi,rcx
  945288:	mov    rdi,r12
  94528b:	call   91d820 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  945290:	mov    rax,QWORD PTR [rsp+0x80]
  945298:	mov    rdx,QWORD PTR [rax+0x38]
  94529c:	mov    rax,QWORD PTR [rax+0x48]
  9452a0:	jmp    944a70 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xe0>
  9452a5:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  9452aa:	mov    BYTE PTR [rsp+0x1b0],r13b
  9452b2:	mov    QWORD PTR [rsp+0x1a0],rax
  9452ba:	lea    rax,[rbp+rbp*2+0x0]
  9452bf:	mov    QWORD PTR [rsp+0x1a8],0x0
  9452cb:	shl    rax,0x3
  9452cf:	movdqa xmm6,XMMWORD PTR [rsp+0x1a0]
  9452d8:	movups XMMWORD PTR [rbx+rax*1-0x360],xmm6
  9452e0:	mov    BYTE PTR [rax+r14*1+0x10],r13b
  9452e5:	movzx  ebp,BYTE PTR [r14+0x180]
  9452ed:	jmp    944a16 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x86>
  9452f2:	movzx  eax,WORD PTR [r13+0xfe2]
  9452fa:	cmp    ax,0xffff
  9452fe:	je     9450ff <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x76f>
  945304:	add    eax,0x1
  945307:	mov    WORD PTR [r13+0xfe2],ax
  94530f:	jmp    9450ff <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x76f>
  945314:	call   920350 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  945319:	mov    rax,QWORD PTR [rsp+0x288]
  945321:	sub    rax,QWORD PTR fs:0x28
  94532a:	je     94533a <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x9aa>
  94532c:	call   e0f90 <__stack_chk_fail@plt>
  945331:	endbr64
  945335:	mov    rbx,rax
  945338:	jmp    945314 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x984>
  94533a:	mov    rdi,rbx
  94533d:	call   e1b80 <_Unwind_Resume@plt>

Disassembly of section .fini:
