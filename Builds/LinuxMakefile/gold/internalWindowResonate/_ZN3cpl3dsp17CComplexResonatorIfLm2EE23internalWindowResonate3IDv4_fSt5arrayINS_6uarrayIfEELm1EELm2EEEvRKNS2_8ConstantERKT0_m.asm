; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

000000000093f7b0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long):
  93f7b0:	endbr64
  93f7b4:	push   r15
  93f7b6:	push   r14
  93f7b8:	push   r13
  93f7ba:	push   r12
  93f7bc:	mov    r12,rdi
  93f7bf:	push   rbp
  93f7c0:	push   rbx
  93f7c1:	sub    rsp,0x298
  93f7c8:	mov    QWORD PTR [rsp+0x80],rsi
  93f7d0:	mov    QWORD PTR [rsp+0x88],rdx
  93f7d8:	mov    QWORD PTR [rsp+0xc8],rcx
  93f7e0:	mov    rax,QWORD PTR fs:0x28
  93f7e9:	mov    QWORD PTR [rsp+0x288],rax
  93f7f1:	xor    eax,eax
  93f7f3:	movzx  r13d,BYTE PTR [rip+0x2df8f7]        # c1f0f2 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  93f7fb:	test   r13b,r13b
  93f7fe:	je     940057 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x8a7>
  93f804:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  93f810:	cmp    QWORD PTR [rax-0x1d8],0x0
  93f818:	mov    rbx,rax
  93f81b:	lea    r14,[rax-0x360]
  93f822:	je     93f84c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x9c>
  93f824:	movzx  ebp,BYTE PTR [r14+0x180]
  93f82c:	cmp    bpl,0xf
  93f830:	jbe    9400c5 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x915>
  93f836:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  93f842:	add    ebp,0x1
  93f845:	mov    BYTE PTR [rax-0x1e0],bpl
  93f84c:	mov    rsi,QWORD PTR [rsp+0x80]
  93f854:	mov    r9,QWORD PTR [r12+0x8]
  93f859:	mov    rdi,QWORD PTR [r12]
  93f85d:	mov    rdx,QWORD PTR [rsi+0x38]
  93f861:	mov    rax,QWORD PTR [rsi+0x48]
  93f865:	mov    rcx,r9
  93f868:	sub    rcx,rdi
  93f86b:	mov    r8,rdx
  93f86e:	sar    rcx,0x2
  93f872:	imul   r8,rax
  93f876:	lea    rsi,[r8*4+0x0]
  93f87e:	cmp    rcx,rsi
  93f881:	jb     9400a5 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x8f5>
  93f887:	cmp    rsi,rcx
  93f88a:	jb     94003d <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x88d>
  93f890:	mov    rcx,QWORD PTR [rsp+0x80]
  93f898:	lea    rsi,[rax+rax*1]
  93f89c:	imul   rdx,rsi
  93f8a0:	cmp    QWORD PTR [rcx+0x40],0x0
  93f8a5:	je     93fedf <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x72f>
  93f8ab:	lea    r8,[rsi+rax*1]
  93f8af:	lea    rcx,[rax*4+0x0]
  93f8b7:	xor    r14d,r14d
  93f8ba:	mov    QWORD PTR [rsp+0x70],0x0
  93f8c3:	lea    rbx,[r8*4+0x0]
  93f8cb:	lea    rdi,[rcx+rdx*1]
  93f8cf:	mov    QWORD PTR [rsp+0x198],rcx
  93f8d7:	mov    r13,rcx
  93f8da:	mov    QWORD PTR [rsp+0x60],rbx
  93f8df:	mov    rbx,rax
  93f8e2:	add    rax,rcx
  93f8e5:	shl    rax,0x2
  93f8e9:	shl    rbx,0x4
  93f8ed:	mov    QWORD PTR [rsp+0x50],rax
  93f8f2:	lea    rax,[rdx*4+0x0]
  93f8fa:	mov    QWORD PTR [rsp+0x68],rax
  93f8ff:	mov    rax,rdi
  93f902:	sub    rax,rsi
  93f905:	mov    QWORD PTR [rsp+0x58],rbx
  93f90a:	shl    rax,0x2
  93f90e:	mov    QWORD PTR [rsp+0x98],rax
  93f916:	lea    rax,[rdi*4+0x0]
  93f91e:	mov    QWORD PTR [rsp+0x90],rax
  93f926:	lea    rax,[rsp+0x1c0]
  93f92e:	mov    QWORD PTR [rsp+0xc0],rax
  93f936:	lea    rax,[rsp+0x1e0]
  93f93e:	mov    QWORD PTR [rsp+0xa0],rax
  93f946:	lea    rax,[rsp+0x200]
  93f94e:	mov    QWORD PTR [rsp+0xa8],rax
  93f956:	lea    rax,[rsp+0x220]
  93f95e:	mov    QWORD PTR [rsp+0xb0],rax
  93f966:	lea    rax,[rsp+0x240]
  93f96e:	mov    QWORD PTR [rsp+0xb8],rax
  93f976:	lea    rax,[rsp+0x260]
  93f97e:	mov    QWORD PTR [rsp+0x190],rax
  93f986:	cs nop WORD PTR [rax+rax*1+0x0]
  93f990:	mov    rax,QWORD PTR [rsp+0x80]
  93f998:	mov    rbx,QWORD PTR [rax]
  93f99b:	lea    rdi,[rbx+r14*1]
  93f99f:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f9a4:	lea    rdi,[rbx+r13*1]
  93f9a8:	movaps XMMWORD PTR [rsp],xmm0
  93f9ac:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f9b1:	mov    rax,QWORD PTR [rsp+0x198]
  93f9b9:	movaps XMMWORD PTR [rsp+0x10],xmm0
  93f9be:	lea    r15,[rax+r13*1]
  93f9c2:	lea    rdi,[rbx+r15*1]
  93f9c6:	mov    QWORD PTR [rsp+0x78],r15
  93f9cb:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f9d0:	mov    rax,QWORD PTR [rsp+0x60]
  93f9d5:	movaps XMMWORD PTR [rsp+0x20],xmm0
  93f9da:	lea    rdi,[rbx+rax*1]
  93f9de:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f9e3:	mov    rsi,QWORD PTR [rsp+0x58]
  93f9e8:	movaps XMMWORD PTR [rsp+0x30],xmm0
  93f9ed:	lea    rdi,[rbx+rsi*1]
  93f9f1:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f9f6:	mov    rcx,QWORD PTR [rsp+0x50]
  93f9fb:	movaps XMMWORD PTR [rsp+0x40],xmm0
  93fa00:	lea    rdi,[rbx+rcx*1]
  93fa04:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fa09:	mov    rbx,QWORD PTR [r12]
  93fa0d:	mov    rdx,QWORD PTR [rsp+0x88]
  93fa15:	movaps XMMWORD PTR [rsp+0x180],xmm0
  93fa1d:	lea    rdi,[rbx+r14*1]
  93fa21:	mov    rbp,QWORD PTR [rdx]
  93fa24:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fa29:	lea    rdi,[rbx+r13*1]
  93fa2d:	movaps XMMWORD PTR [rsp+0x1c0],xmm0
  93fa35:	movaps XMMWORD PTR [rsp+0x170],xmm0
  93fa3d:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fa42:	lea    rdi,[rbx+r15*1]
  93fa46:	movaps XMMWORD PTR [rsp+0x1e0],xmm0
  93fa4e:	movaps XMMWORD PTR [rsp+0x160],xmm0
  93fa56:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fa5b:	mov    rax,QWORD PTR [rsp+0x60]
  93fa60:	movaps XMMWORD PTR [rsp+0x200],xmm0
  93fa68:	lea    rdi,[rbx+rax*1]
  93fa6c:	movaps XMMWORD PTR [rsp+0x150],xmm0
  93fa74:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fa79:	mov    rsi,QWORD PTR [rsp+0x58]
  93fa7e:	movaps XMMWORD PTR [rsp+0x220],xmm0
  93fa86:	lea    rdi,[rbx+rsi*1]
  93fa8a:	movaps XMMWORD PTR [rsp+0x140],xmm0
  93fa92:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fa97:	mov    rcx,QWORD PTR [rsp+0x50]
  93fa9c:	movaps XMMWORD PTR [rsp+0x240],xmm0
  93faa4:	lea    rdi,[rbx+rcx*1]
  93faa8:	movaps XMMWORD PTR [rsp+0x130],xmm0
  93fab0:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fab5:	mov    rax,QWORD PTR [rsp+0x68]
  93faba:	mov    rdx,QWORD PTR [rsp+0x88]
  93fac2:	movaps XMMWORD PTR [rsp+0x260],xmm0
  93faca:	lea    rdi,[rax+r14*1]
  93face:	mov    r15,QWORD PTR [rdx+0x10]
  93fad2:	movaps XMMWORD PTR [rsp+0x120],xmm0
  93fada:	add    rdi,rbx
  93fadd:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fae2:	mov    rax,QWORD PTR [rsp+0x68]
  93fae7:	movaps XMMWORD PTR [rsp+0x1d0],xmm0
  93faef:	lea    rdi,[rax+r13*1]
  93faf3:	movaps XMMWORD PTR [rsp+0x110],xmm0
  93fafb:	add    rdi,rbx
  93fafe:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fb03:	movaps XMMWORD PTR [rsp+0x1f0],xmm0
  93fb0b:	movaps XMMWORD PTR [rsp+0x100],xmm0
  93fb13:	mov    rax,QWORD PTR [rsp+0x98]
  93fb1b:	lea    rdi,[rax+r14*1]
  93fb1f:	add    rdi,rbx
  93fb22:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fb27:	mov    rax,QWORD PTR [rsp+0x98]
  93fb2f:	movaps XMMWORD PTR [rsp+0x210],xmm0
  93fb37:	lea    rdi,[rax+r13*1]
  93fb3b:	movaps XMMWORD PTR [rsp+0xf0],xmm0
  93fb43:	add    rdi,rbx
  93fb46:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fb4b:	mov    rax,QWORD PTR [rsp+0x90]
  93fb53:	movaps XMMWORD PTR [rsp+0x230],xmm0
  93fb5b:	lea    rdi,[rax+r14*1]
  93fb5f:	movaps XMMWORD PTR [rsp+0xe0],xmm0
  93fb67:	add    rdi,rbx
  93fb6a:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fb6f:	mov    rax,QWORD PTR [rsp+0x90]
  93fb77:	movaps XMMWORD PTR [rsp+0x250],xmm0
  93fb7f:	lea    rdi,[rax+r13*1]
  93fb83:	movaps XMMWORD PTR [rsp+0xd0],xmm0
  93fb8b:	add    rdi,rbx
  93fb8e:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93fb93:	xor    eax,eax
  93fb95:	cmp    QWORD PTR [rsp+0xc8],0x0
  93fb9e:	movaps xmm1,XMMWORD PTR [rsp+0xd0]
  93fba6:	movaps xmm4,XMMWORD PTR [rsp+0xe0]
  93fbae:	movaps xmm3,XMMWORD PTR [rsp+0xf0]
  93fbb6:	movaps xmm2,xmm0
  93fbb9:	movaps XMMWORD PTR [rsp+0x270],xmm0
  93fbc1:	movaps xmm6,XMMWORD PTR [rsp+0x100]
  93fbc9:	movaps xmm5,XMMWORD PTR [rsp+0x110]
  93fbd1:	movaps xmm8,XMMWORD PTR [rsp+0x120]
  93fbda:	movaps xmm7,XMMWORD PTR [rsp+0x130]
  93fbe2:	movaps xmm10,XMMWORD PTR [rsp+0x140]
  93fbeb:	movaps xmm9,XMMWORD PTR [rsp+0x150]
  93fbf4:	movaps xmm12,XMMWORD PTR [rsp+0x160]
  93fbfd:	movaps xmm11,XMMWORD PTR [rsp+0x170]
  93fc06:	je     93fdb0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x600>
  93fc0c:	movaps xmm13,XMMWORD PTR [rsp+0x180]
  93fc15:	mov    rdx,QWORD PTR [rsp+0xc8]
  93fc1d:	nop    DWORD PTR [rax]
  93fc20:	movss  xmm0,DWORD PTR [rbp+rax*4+0x0]
  93fc26:	call   9205e0 <float __vector(4) cpl::simd::broadcast<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*) [clone .isra.0]>
  93fc2b:	movaps xmm14,XMMWORD PTR [rsp]
  93fc30:	movaps xmm15,XMMWORD PTR [rsp+0x10]
  93fc36:	mulps  xmm14,xmm11
  93fc3a:	mulps  xmm15,xmm12
  93fc3e:	mulps  xmm11,XMMWORD PTR [rsp+0x10]
  93fc44:	mulps  xmm12,XMMWORD PTR [rsp]
  93fc49:	subps  xmm14,xmm15
  93fc4d:	movaps xmm15,XMMWORD PTR [rsp+0x30]
  93fc53:	addps  xmm12,xmm11
  93fc57:	mulps  xmm15,xmm10
  93fc5b:	mulps  xmm10,XMMWORD PTR [rsp+0x20]
  93fc61:	movaps xmm11,xmm14
  93fc65:	movaps xmm14,XMMWORD PTR [rsp+0x20]
  93fc6b:	addps  xmm11,xmm0
  93fc6f:	mulps  xmm14,xmm9
  93fc73:	mulps  xmm9,XMMWORD PTR [rsp+0x30]
  93fc79:	subps  xmm14,xmm15
  93fc7d:	movaps xmm15,xmm8
  93fc81:	addps  xmm10,xmm9
  93fc85:	mulps  xmm15,xmm13
  93fc89:	mulps  xmm8,XMMWORD PTR [rsp+0x40]
  93fc8f:	movaps xmm9,xmm14
  93fc93:	movaps xmm14,XMMWORD PTR [rsp+0x40]
  93fc99:	addps  xmm9,xmm0
  93fc9d:	mulps  xmm14,xmm7
  93fca1:	mulps  xmm7,xmm13
  93fca5:	subps  xmm14,xmm15
  93fca9:	addps  xmm8,xmm7
  93fcad:	movaps xmm7,xmm0
  93fcb0:	movss  xmm0,DWORD PTR [r15+rax*4]
  93fcb6:	add    rax,0x1
  93fcba:	call   9205e0 <float __vector(4) cpl::simd::broadcast<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*) [clone .isra.0]>
  93fcbf:	addps  xmm7,xmm14
  93fcc3:	movaps xmm15,XMMWORD PTR [rsp+0x10]
  93fcc9:	movaps xmm14,XMMWORD PTR [rsp]
  93fcce:	mulps  xmm15,xmm6
  93fcd2:	mulps  xmm14,xmm5
  93fcd6:	mulps  xmm6,XMMWORD PTR [rsp]
  93fcda:	mulps  xmm5,XMMWORD PTR [rsp+0x10]
  93fcdf:	subps  xmm14,xmm15
  93fce3:	movaps xmm15,XMMWORD PTR [rsp+0x30]
  93fce9:	addps  xmm6,xmm5
  93fcec:	mulps  xmm15,xmm4
  93fcf0:	mulps  xmm4,XMMWORD PTR [rsp+0x20]
  93fcf5:	movaps xmm5,xmm14
  93fcf9:	movaps xmm14,XMMWORD PTR [rsp+0x20]
  93fcff:	addps  xmm5,xmm0
  93fd02:	mulps  xmm14,xmm3
  93fd06:	mulps  xmm3,XMMWORD PTR [rsp+0x30]
  93fd0b:	subps  xmm14,xmm15
  93fd0f:	movaps xmm15,xmm13
  93fd13:	addps  xmm4,xmm3
  93fd16:	mulps  xmm15,xmm2
  93fd1a:	mulps  xmm2,XMMWORD PTR [rsp+0x40]
  93fd1f:	movaps xmm3,xmm14
  93fd23:	movaps xmm14,XMMWORD PTR [rsp+0x40]
  93fd29:	addps  xmm3,xmm0
  93fd2c:	mulps  xmm14,xmm1
  93fd30:	mulps  xmm1,xmm13
  93fd34:	subps  xmm14,xmm15
  93fd38:	addps  xmm2,xmm1
  93fd3b:	movaps xmm1,xmm0
  93fd3e:	addps  xmm1,xmm14
  93fd42:	cmp    rdx,rax
  93fd45:	jne    93fc20 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x470>
  93fd4b:	movaps XMMWORD PTR [rsp+0x1e0],xmm12
  93fd54:	movaps XMMWORD PTR [rsp+0x1c0],xmm11
  93fd5d:	movaps XMMWORD PTR [rsp+0x220],xmm10
  93fd66:	movaps XMMWORD PTR [rsp+0x200],xmm9
  93fd6f:	movaps XMMWORD PTR [rsp+0x260],xmm8
  93fd78:	movaps XMMWORD PTR [rsp+0x240],xmm7
  93fd80:	movaps XMMWORD PTR [rsp+0x1f0],xmm6
  93fd88:	movaps XMMWORD PTR [rsp+0x1d0],xmm5
  93fd90:	movaps XMMWORD PTR [rsp+0x230],xmm4
  93fd98:	movaps XMMWORD PTR [rsp+0x210],xmm3
  93fda0:	movaps XMMWORD PTR [rsp+0x270],xmm2
  93fda8:	movaps XMMWORD PTR [rsp+0x250],xmm1
  93fdb0:	xor    r15d,r15d
  93fdb3:	mov    rdi,rbx
  93fdb6:	xor    ebp,ebp
  93fdb8:	mov    rbx,r15
  93fdbb:	mov    r15,QWORD PTR [rsp+0x190]
  93fdc3:	mov    rax,QWORD PTR [rsp+0xc0]
  93fdcb:	lea    rdx,[r14+rbp*1]
  93fdcf:	add    rdi,rdx
  93fdd2:	movaps xmm2,XMMWORD PTR [rax+rbx*1]
  93fdd6:	movaps xmm0,xmm2
  93fdd9:	movaps XMMWORD PTR [rsp],xmm2
  93fddd:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  93fde2:	mov    rax,QWORD PTR [rsp+0xa0]
  93fdea:	lea    rdi,[r13+rbp*1+0x0]
  93fdef:	add    rdi,QWORD PTR [r12]
  93fdf3:	movaps xmm4,XMMWORD PTR [rax+rbx*1]
  93fdf7:	movaps xmm0,xmm4
  93fdfa:	movaps XMMWORD PTR [rsp],xmm4
  93fdfe:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  93fe03:	mov    rax,QWORD PTR [rsp+0xa8]
  93fe0b:	movaps xmm6,XMMWORD PTR [rax+rbx*1]
  93fe0f:	mov    rax,QWORD PTR [rsp+0x78]
  93fe14:	movaps xmm0,xmm6
  93fe17:	lea    rdi,[rax+rbp*1]
  93fe1b:	add    rdi,QWORD PTR [r12]
  93fe1f:	movaps XMMWORD PTR [rsp],xmm6
  93fe23:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  93fe28:	mov    rax,QWORD PTR [rsp+0xb0]
  93fe30:	movaps xmm2,XMMWORD PTR [rax+rbx*1]
  93fe34:	mov    rax,QWORD PTR [rsp+0x60]
  93fe39:	movaps xmm0,xmm2
  93fe3c:	lea    rdi,[rax+rbp*1]
  93fe40:	add    rdi,QWORD PTR [r12]
  93fe44:	movaps XMMWORD PTR [rsp],xmm2
  93fe48:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  93fe4d:	mov    rax,QWORD PTR [rsp+0xb8]
  93fe55:	movaps xmm4,XMMWORD PTR [rax+rbx*1]
  93fe59:	mov    rax,QWORD PTR [rsp+0x58]
  93fe5e:	movaps xmm0,xmm4
  93fe61:	lea    rdi,[rax+rbp*1]
  93fe65:	add    rdi,QWORD PTR [r12]
  93fe69:	movaps XMMWORD PTR [rsp],xmm4
  93fe6d:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  93fe72:	mov    rax,QWORD PTR [rsp+0x50]
  93fe77:	movaps xmm0,XMMWORD PTR [r15+rbx*1]
  93fe7c:	add    rbx,0x10
  93fe80:	lea    rdi,[rax+rbp*1]
  93fe84:	add    rdi,QWORD PTR [r12]
  93fe88:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  93fe8d:	cmp    rbx,0x20
  93fe91:	je     93fea8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x6f8>
  93fe93:	mov    rax,QWORD PTR [rsp+0x68]
  93fe98:	mov    rdi,QWORD PTR [r12]
  93fe9c:	add    rbp,rax
  93fe9f:	jmp    93fdc3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x613>
  93fea4:	nop    DWORD PTR [rax+0x0]
  93fea8:	add    QWORD PTR [rsp+0x70],0x4
  93feae:	add    r14,0x10
  93feb2:	mov    rax,QWORD PTR [rsp+0x70]
  93feb7:	add    r13,0x10
  93febb:	mov    rcx,QWORD PTR [rsp+0x80]
  93fec3:	add    QWORD PTR [rsp+0x60],0x10
  93fec9:	add    QWORD PTR [rsp+0x58],0x10
  93fecf:	add    QWORD PTR [rsp+0x50],0x10
  93fed5:	cmp    rax,QWORD PTR [rcx+0x40]
  93fed9:	jb     93f990 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x1e0>
  93fedf:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  93feeb:	lea    rbp,[rax-0x360]
  93fef2:	mov    r12,rax
  93fef5:	mov    r13,QWORD PTR [rbp+0x188]
  93fefc:	test   r13,r13
  93feff:	je     93ff1f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x76f>
  93ff01:	movzx  eax,BYTE PTR [rbp+0x180]
  93ff08:	lea    ebx,[rax-0x1]
  93ff0b:	mov    BYTE PTR [rbp+0x180],bl
  93ff11:	cmp    bl,0xf
  93ff14:	jbe    93ff48 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x798>
  93ff16:	add    WORD PTR [r13+0xfe2],0x1
  93ff1f:	mov    rax,QWORD PTR [rsp+0x288]
  93ff27:	sub    rax,QWORD PTR fs:0x28
  93ff30:	jne    94014c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x99c>
  93ff36:	add    rsp,0x298
  93ff3d:	pop    rbx
  93ff3e:	pop    rbp
  93ff3f:	pop    r12
  93ff41:	pop    r13
  93ff43:	pop    r14
  93ff45:	pop    r15
  93ff47:	ret
  93ff48:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  93ff4d:	movzx  r15d,bl
  93ff51:	lea    rdx,[r15+r15*2]
  93ff55:	mov    rcx,rax
  93ff58:	movq   xmm0,rax
  93ff5d:	shl    rdx,0x3
  93ff61:	movdqu xmm1,XMMWORD PTR [r12+rdx*1-0x360]
  93ff6b:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  93ff70:	movq   xmm6,rcx
  93ff75:	movzx  r12d,BYTE PTR [rbp+0x181]
  93ff7d:	punpcklqdq xmm0,xmm6
  93ff81:	psubq  xmm0,xmm1
  93ff85:	cmp    r12b,bl
  93ff88:	jae    93ffa8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x7f8>
  93ff8a:	lea    eax,[r15-0x1]
  93ff8e:	movdqa xmm1,xmm0
  93ff92:	cdqe
  93ff94:	lea    rax,[rax+rax*2]
  93ff98:	movq   xmm2,QWORD PTR [rbp+rax*8+0x8]
  93ff9e:	paddq  xmm1,xmm2
  93ffa2:	movq   QWORD PTR [rbp+rax*8+0x8],xmm1
  93ffa8:	movzx  ebp,WORD PTR [r13+0xfe0]
  93ffb0:	cmp    bp,0x7f
  93ffb4:	je     940112 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x962>
  93ffba:	movaps XMMWORD PTR [rsp],xmm0
  93ffbe:	sub    ebx,r12d
  93ffc1:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  93ffcd:	lea    rdx,[r15+r15*2]
  93ffd1:	movdqa xmm0,XMMWORD PTR [rsp]
  93ffd6:	mov    DWORD PTR [rsp+0x1b8],0x0
  93ffe1:	movups XMMWORD PTR [rsp+0x1a8],xmm0
  93ffe9:	lea    rdx,[rax+rdx*8-0x360]
  93fff1:	mov    rax,QWORD PTR [rdx]
  93fff4:	mov    QWORD PTR [rsp+0x1a0],rax
  93fffc:	movzx  eax,BYTE PTR [rdx+0x10]
  940000:	lea    edx,[rbp+0x1]
  940003:	movdqa xmm6,XMMWORD PTR [rsp+0x1a0]
  94000c:	mov    WORD PTR [r13+0xfe0],dx
  940014:	mov    ah,bl
  940016:	mov    WORD PTR [rsp+0x1bc],ax
  94001e:	movzx  eax,bp
  940021:	shl    rax,0x5
  940025:	add    rax,r13
  940028:	movups XMMWORD PTR [rax],xmm6
  94002b:	movdqu xmm6,XMMWORD PTR [rsp+0x1ae]
  940034:	movups XMMWORD PTR [rax+0xe],xmm6
  940038:	jmp    93ff1f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x76f>
  94003d:	shl    r8,0x4
  940041:	add    rdi,r8
  940044:	cmp    r9,rdi
  940047:	je     93f890 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xe0>
  94004d:	mov    QWORD PTR [r12+0x8],rdi
  940052:	jmp    93f890 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xe0>
  940057:	mov    eax,0x1
  94005c:	lock xadd DWORD PTR [rip+0x2d3afc],eax        # c13b60 <cpl::Profiling::registerRegion(char const*)::counter>
  940064:	mov    edx,0x1
  940069:	mov    r13d,0x1
  94006f:	add    eax,0x2
  940072:	cmp    eax,0xfe
  940077:	ja     940092 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x8e2>
  940079:	lea    rdx,[rip+0x2d1fa0]        # c12020 <cpl::Profiling::regions>
  940080:	mov    ecx,eax
  940082:	lea    rsi,[rip+0xa80bd]        # 9e8146 <_IO_stdin_used+0x18146>
  940089:	mov    r13d,eax
  94008c:	mov    QWORD PTR [rdx+rcx*8],rsi
  940090:	mov    edx,eax
  940092:	xor    eax,eax
  940094:	lock cmpxchg BYTE PTR [rip+0x2df056],dl        # c1f0f2 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  94009c:	cmovne r13d,eax
  9400a0:	jmp    93f804 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x54>
  9400a5:	sub    rsi,rcx
  9400a8:	mov    rdi,r12
  9400ab:	call   91d820 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  9400b0:	mov    rax,QWORD PTR [rsp+0x80]
  9400b8:	mov    rdx,QWORD PTR [rax+0x38]
  9400bc:	mov    rax,QWORD PTR [rax+0x48]
  9400c0:	jmp    93f890 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xe0>
  9400c5:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  9400ca:	mov    BYTE PTR [rsp+0x1b0],r13b
  9400d2:	mov    QWORD PTR [rsp+0x1a0],rax
  9400da:	lea    rax,[rbp+rbp*2+0x0]
  9400df:	mov    QWORD PTR [rsp+0x1a8],0x0
  9400eb:	shl    rax,0x3
  9400ef:	movdqa xmm6,XMMWORD PTR [rsp+0x1a0]
  9400f8:	movups XMMWORD PTR [rbx+rax*1-0x360],xmm6
  940100:	mov    BYTE PTR [rax+r14*1+0x10],r13b
  940105:	movzx  ebp,BYTE PTR [r14+0x180]
  94010d:	jmp    93f836 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x86>
  940112:	movzx  eax,WORD PTR [r13+0xfe2]
  94011a:	cmp    ax,0xffff
  94011e:	je     93ff1f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x76f>
  940124:	add    eax,0x1
  940127:	mov    WORD PTR [r13+0xfe2],ax
  94012f:	jmp    93ff1f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x76f>
  940134:	call   920350 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  940139:	mov    rax,QWORD PTR [rsp+0x288]
  940141:	sub    rax,QWORD PTR fs:0x28
  94014a:	je     94015a <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x9aa>
  94014c:	call   e0f90 <__stack_chk_fail@plt>
  940151:	endbr64
  940155:	mov    rbx,rax
  940158:	jmp    940134 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x984>
  94015a:	mov    rdi,rbx
  94015d:	call   e1b80 <_Unwind_Resume@plt>

Disassembly of section .fini:
