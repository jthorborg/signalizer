; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

000000000094d7e0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long):
  94d7e0:	endbr64
  94d7e4:	push   r15
  94d7e6:	push   r14
  94d7e8:	push   r13
  94d7ea:	push   r12
  94d7ec:	mov    r12,rdi
  94d7ef:	push   rbp
  94d7f0:	mov    rbp,rcx
  94d7f3:	push   rbx
  94d7f4:	mov    rbx,rdx
  94d7f7:	sub    rsp,0xa8
  94d7fe:	mov    QWORD PTR [rsp],rsi
  94d802:	movzx  r15d,BYTE PTR [rip+0x2d18ce]        # c1f0d8 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  94d80a:	test   r15b,r15b
  94d80d:	je     94dcb4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4d4>
  94d813:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94d81f:	cmp    QWORD PTR [rax-0x1d8],0x0
  94d827:	mov    r14,rax
  94d82a:	lea    rdx,[rax-0x360]
  94d831:	je     94d85c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x7c>
  94d833:	movzx  r13d,BYTE PTR [rax-0x1e0]
  94d83b:	cmp    r13b,0xf
  94d83f:	jbe    94dd1a <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x53a>
  94d845:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94d851:	add    r13d,0x1
  94d855:	mov    BYTE PTR [rax-0x1e0],r13b
  94d85c:	mov    rdx,QWORD PTR [rsp]
  94d860:	mov    r8,QWORD PTR [r12+0x8]
  94d865:	mov    rcx,QWORD PTR [r12]
  94d869:	mov    rax,QWORD PTR [rdx+0x48]
  94d86d:	mov    rdi,QWORD PTR [rdx+0x38]
  94d871:	mov    rdx,r8
  94d874:	sub    rdx,rcx
  94d877:	imul   rdi,rax
  94d87b:	sar    rdx,0x2
  94d87f:	lea    rsi,[rdi*4+0x0]
  94d887:	cmp    rdx,rsi
  94d88a:	jb     94dd02 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x522>
  94d890:	cmp    rsi,rdx
  94d893:	jb     94dc9a <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4ba>
  94d899:	mov    rdx,QWORD PTR [rsp]
  94d89d:	mov    r11,QWORD PTR [rdx+0x40]
  94d8a1:	test   r11,r11
  94d8a4:	je     94db54 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x374>
  94d8aa:	mov    r15,QWORD PTR [r12]
  94d8ae:	mov    rsi,QWORD PTR [rdx]
  94d8b1:	xor    r12d,r12d
  94d8b4:	lea    rdx,[rax*8+0x0]
  94d8bc:	lea    r13,[rax*4+0x0]
  94d8c4:	shl    rax,0x4
  94d8c8:	mov    QWORD PTR [rsp+0x68],r11
  94d8cd:	mov    rbx,QWORD PTR [rbx]
  94d8d0:	lea    rcx,[rsi+rdx*1]
  94d8d4:	add    rax,r15
  94d8d7:	lea    r14,[r15+rdx*1]
  94d8db:	mov    QWORD PTR [rsp+0x70],rsi
  94d8e0:	mov    QWORD PTR [rsp+0x58],rax
  94d8e5:	lea    rdi,[rcx+rdx*1]
  94d8e9:	lea    rax,[r14+rdx*1]
  94d8ed:	mov    QWORD PTR [rsp+0x50],rdi
  94d8f2:	mov    QWORD PTR [rsp+0x60],rax
  94d8f7:	mov    QWORD PTR [rsp],0x0
  94d8ff:	mov    QWORD PTR [rsp+0x78],rbp
  94d904:	mov    rbp,r13
  94d907:	mov    r13,rcx
  94d90a:	nop    WORD PTR [rax+rax*1+0x0]
  94d910:	mov    rax,QWORD PTR [rsp+0x70]
  94d915:	lea    rdi,[r15+r12*1]
  94d919:	movss  xmm0,DWORD PTR [rax+r12*1]
  94d91f:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94d924:	movaps xmm8,xmm0
  94d928:	movss  xmm0,DWORD PTR [rax+rbp*1]
  94d92d:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94d932:	movaps xmm9,xmm0
  94d936:	movss  xmm0,DWORD PTR [r13+r12*1+0x0]
  94d93d:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94d942:	movaps xmm10,xmm0
  94d946:	movss  xmm0,DWORD PTR [r13+rbp*1+0x0]
  94d94d:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94d952:	mov    rax,QWORD PTR [rsp+0x50]
  94d957:	movaps xmm11,xmm0
  94d95b:	movss  xmm0,DWORD PTR [rax+r12*1]
  94d961:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94d966:	movaps xmm12,xmm0
  94d96a:	movss  xmm0,DWORD PTR [rax+rbp*1]
  94d96f:	lea    rax,[r15+rbp*1]
  94d973:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94d978:	movaps xmm13,xmm0
  94d97c:	movss  xmm0,DWORD PTR [r15+r12*1]
  94d982:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94d987:	mov    QWORD PTR [rsp+0x30],rax
  94d98c:	lea    rax,[r14+r12*1]
  94d990:	movaps xmm6,xmm0
  94d993:	movss  xmm0,DWORD PTR [r15+rbp*1]
  94d999:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94d99e:	mov    QWORD PTR [rsp+0x28],rax
  94d9a3:	lea    rax,[r14+rbp*1]
  94d9a7:	movaps xmm5,xmm0
  94d9aa:	movss  xmm0,DWORD PTR [r14+r12*1]
  94d9b0:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94d9b5:	mov    QWORD PTR [rsp+0x10],rax
  94d9ba:	movaps xmm4,xmm0
  94d9bd:	movss  xmm0,DWORD PTR [r14+rbp*1]
  94d9c3:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94d9c8:	mov    rcx,QWORD PTR [rsp+0x60]
  94d9cd:	mov    rax,QWORD PTR [rsp+0x58]
  94d9d2:	movaps xmm3,xmm0
  94d9d5:	movss  xmm0,DWORD PTR [rcx+r12*1]
  94d9db:	lea    rdx,[rax+r12*1]
  94d9df:	add    rax,rbp
  94d9e2:	mov    QWORD PTR [rsp+0x18],rdx
  94d9e7:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94d9ec:	mov    QWORD PTR [rsp+0x20],rax
  94d9f1:	movaps xmm2,xmm0
  94d9f4:	movss  xmm0,DWORD PTR [rcx+rbp*1]
  94d9f9:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94d9fe:	mov    rsi,QWORD PTR [rsp+0x78]
  94da03:	movaps xmm1,xmm0
  94da06:	test   rsi,rsi
  94da09:	je     94dab0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x2d0>
  94da0f:	xor    eax,eax
  94da11:	nop    DWORD PTR [rax+0x0]
  94da18:	movaps xmm7,xmm6
  94da1b:	movaps xmm14,xmm5
  94da1f:	movss  xmm0,DWORD PTR [rbx+rax*4]
  94da24:	add    rax,0x1
  94da28:	mulss  xmm14,xmm9
  94da2d:	call   9205f0 <float cpl::simd::broadcast<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94da32:	mulss  xmm7,xmm8
  94da37:	mulss  xmm6,xmm9
  94da3c:	mulss  xmm5,xmm8
  94da41:	subss  xmm7,xmm14
  94da46:	movaps xmm14,xmm11
  94da4a:	mulss  xmm14,xmm3
  94da4f:	addss  xmm5,xmm6
  94da53:	mulss  xmm3,xmm10
  94da58:	movaps xmm6,xmm7
  94da5b:	movaps xmm7,xmm4
  94da5e:	mulss  xmm7,xmm10
  94da63:	addss  xmm6,xmm0
  94da67:	mulss  xmm4,xmm11
  94da6c:	subss  xmm7,xmm14
  94da71:	movaps xmm14,xmm1
  94da75:	addss  xmm3,xmm4
  94da79:	mulss  xmm14,xmm13
  94da7e:	mulss  xmm1,xmm12
  94da83:	movaps xmm4,xmm7
  94da86:	movaps xmm7,xmm2
  94da89:	mulss  xmm7,xmm12
  94da8e:	addss  xmm4,xmm0
  94da92:	mulss  xmm2,xmm13
  94da97:	subss  xmm7,xmm14
  94da9c:	addss  xmm1,xmm2
  94daa0:	movaps xmm2,xmm0
  94daa3:	addss  xmm2,xmm7
  94daa7:	cmp    rsi,rax
  94daaa:	jne    94da18 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x238>
  94dab0:	movaps xmm0,xmm6
  94dab3:	movss  DWORD PTR [rsp+0x3c],xmm1
  94dab9:	add    r12,0x4
  94dabd:	add    rbp,0x4
  94dac1:	movss  DWORD PTR [rsp+0x40],xmm2
  94dac7:	movss  DWORD PTR [rsp+0x48],xmm3
  94dacd:	movss  DWORD PTR [rsp+0x4c],xmm4
  94dad3:	movss  DWORD PTR [rsp+0x44],xmm5
  94dad9:	call   85fb30 <cpl::simd::store(float*, float)>
  94dade:	movss  xmm5,DWORD PTR [rsp+0x44]
  94dae4:	mov    rdi,QWORD PTR [rsp+0x30]
  94dae9:	movaps xmm0,xmm5
  94daec:	call   85fb30 <cpl::simd::store(float*, float)>
  94daf1:	movss  xmm4,DWORD PTR [rsp+0x4c]
  94daf7:	mov    rdi,QWORD PTR [rsp+0x28]
  94dafc:	movaps xmm0,xmm4
  94daff:	call   85fb30 <cpl::simd::store(float*, float)>
  94db04:	movss  xmm3,DWORD PTR [rsp+0x48]
  94db0a:	mov    rdi,QWORD PTR [rsp+0x10]
  94db0f:	movaps xmm0,xmm3
  94db12:	call   85fb30 <cpl::simd::store(float*, float)>
  94db17:	movss  xmm2,DWORD PTR [rsp+0x40]
  94db1d:	mov    rdi,QWORD PTR [rsp+0x18]
  94db22:	movaps xmm0,xmm2
  94db25:	call   85fb30 <cpl::simd::store(float*, float)>
  94db2a:	movss  xmm1,DWORD PTR [rsp+0x3c]
  94db30:	mov    rdi,QWORD PTR [rsp+0x20]
  94db35:	movaps xmm0,xmm1
  94db38:	call   85fb30 <cpl::simd::store(float*, float)>
  94db3d:	add    QWORD PTR [rsp],0x1
  94db42:	mov    rdx,QWORD PTR [rsp+0x68]
  94db47:	mov    rax,QWORD PTR [rsp]
  94db4b:	cmp    rax,rdx
  94db4e:	jne    94d910 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x130>
  94db54:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94db60:	lea    rbp,[rax-0x360]
  94db67:	mov    r12,rax
  94db6a:	mov    r13,QWORD PTR [rbp+0x188]
  94db71:	test   r13,r13
  94db74:	je     94db94 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3b4>
  94db76:	movzx  eax,BYTE PTR [rbp+0x180]
  94db7d:	lea    ebx,[rax-0x1]
  94db80:	mov    BYTE PTR [rbp+0x180],bl
  94db86:	cmp    bl,0xf
  94db89:	jbe    94dba6 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3c6>
  94db8b:	add    WORD PTR [r13+0xfe2],0x1
  94db94:	add    rsp,0xa8
  94db9b:	pop    rbx
  94db9c:	pop    rbp
  94db9d:	pop    r12
  94db9f:	pop    r13
  94dba1:	pop    r14
  94dba3:	pop    r15
  94dba5:	ret
  94dba6:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  94dbab:	movzx  r14d,bl
  94dbaf:	lea    rdx,[r14+r14*2]
  94dbb3:	mov    rcx,rax
  94dbb6:	movq   xmm0,rax
  94dbbb:	shl    rdx,0x3
  94dbbf:	movdqu xmm1,XMMWORD PTR [r12+rdx*1-0x360]
  94dbc9:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  94dbce:	movq   xmm3,rcx
  94dbd3:	movzx  r12d,BYTE PTR [rbp+0x181]
  94dbdb:	punpcklqdq xmm0,xmm3
  94dbdf:	psubq  xmm0,xmm1
  94dbe3:	cmp    r12b,bl
  94dbe6:	jae    94dc06 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x426>
  94dbe8:	lea    eax,[r14-0x1]
  94dbec:	movdqa xmm1,xmm0
  94dbf0:	cdqe
  94dbf2:	lea    rax,[rax+rax*2]
  94dbf6:	movq   xmm2,QWORD PTR [rbp+rax*8+0x8]
  94dbfc:	paddq  xmm1,xmm2
  94dc00:	movq   QWORD PTR [rbp+rax*8+0x8],xmm1
  94dc06:	movzx  ebp,WORD PTR [r13+0xfe0]
  94dc0e:	cmp    bp,0x7f
  94dc12:	je     94dd72 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x592>
  94dc18:	movaps XMMWORD PTR [rsp],xmm0
  94dc1c:	sub    ebx,r12d
  94dc1f:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94dc2b:	lea    rdx,[r14+r14*2]
  94dc2f:	movdqa xmm0,XMMWORD PTR [rsp]
  94dc34:	mov    DWORD PTR [rsp+0x98],0x0
  94dc3f:	movups XMMWORD PTR [rsp+0x88],xmm0
  94dc47:	lea    rdx,[rax+rdx*8-0x360]
  94dc4f:	mov    rax,QWORD PTR [rdx]
  94dc52:	mov    QWORD PTR [rsp+0x80],rax
  94dc5a:	movzx  eax,BYTE PTR [rdx+0x10]
  94dc5e:	lea    edx,[rbp+0x1]
  94dc61:	shl    rbp,0x5
  94dc65:	movdqa xmm5,XMMWORD PTR [rsp+0x80]
  94dc6e:	mov    WORD PTR [r13+0xfe0],dx
  94dc76:	mov    ah,bl
  94dc78:	mov    WORD PTR [rsp+0x9c],ax
  94dc80:	lea    rax,[r13+rbp*1+0x0]
  94dc85:	movdqu xmm3,XMMWORD PTR [rsp+0x8e]
  94dc8e:	movups XMMWORD PTR [rax],xmm5
  94dc91:	movups XMMWORD PTR [rax+0xe],xmm3
  94dc95:	jmp    94db94 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3b4>
  94dc9a:	shl    rdi,0x4
  94dc9e:	add    rcx,rdi
  94dca1:	cmp    r8,rcx
  94dca4:	je     94d899 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb9>
  94dcaa:	mov    QWORD PTR [r12+0x8],rcx
  94dcaf:	jmp    94d899 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb9>
  94dcb4:	mov    eax,0x1
  94dcb9:	lock xadd DWORD PTR [rip+0x2c5e9f],eax        # c13b60 <cpl::Profiling::registerRegion(char const*)::counter>
  94dcc1:	mov    edx,0x1
  94dcc6:	mov    r15d,0x1
  94dccc:	add    eax,0x2
  94dccf:	cmp    eax,0xfe
  94dcd4:	ja     94dcef <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x50f>
  94dcd6:	lea    rdx,[rip+0x2c4343]        # c12020 <cpl::Profiling::regions>
  94dcdd:	mov    ecx,eax
  94dcdf:	lea    rdi,[rip+0x9a460]        # 9e8146 <_IO_stdin_used+0x18146>
  94dce6:	mov    r15d,eax
  94dce9:	mov    QWORD PTR [rdx+rcx*8],rdi
  94dced:	mov    edx,eax
  94dcef:	xor    eax,eax
  94dcf1:	lock cmpxchg BYTE PTR [rip+0x2d13df],dl        # c1f0d8 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  94dcf9:	cmovne r15d,eax
  94dcfd:	jmp    94d813 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x33>
  94dd02:	sub    rsi,rdx
  94dd05:	mov    rdi,r12
  94dd08:	call   91d820 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  94dd0d:	mov    rax,QWORD PTR [rsp]
  94dd11:	mov    rax,QWORD PTR [rax+0x48]
  94dd15:	jmp    94d899 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb9>
  94dd1a:	mov    QWORD PTR [rsp+0x10],rdx
  94dd1f:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  94dd24:	mov    rdx,QWORD PTR [rsp+0x10]
  94dd29:	mov    QWORD PTR [rsp+0x88],0x0
  94dd35:	mov    QWORD PTR [rsp+0x80],rax
  94dd3d:	lea    rax,[r13+r13*2+0x0]
  94dd42:	movdqa xmm5,XMMWORD PTR [rsp+0x80]
  94dd4b:	shl    rax,0x3
  94dd4f:	mov    BYTE PTR [rsp+0x90],r15b
  94dd57:	movups XMMWORD PTR [r14+rax*1-0x360],xmm5
  94dd60:	mov    BYTE PTR [rax+rdx*1+0x10],r15b
  94dd65:	movzx  r13d,BYTE PTR [rdx+0x180]
  94dd6d:	jmp    94d845 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x65>
  94dd72:	movzx  eax,WORD PTR [r13+0xfe2]
  94dd7a:	cmp    ax,0xffff
  94dd7e:	je     94db94 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3b4>
  94dd84:	add    eax,0x1
  94dd87:	mov    WORD PTR [r13+0xfe2],ax
  94dd8f:	jmp    94db94 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3b4>
  94dd94:	endbr64
  94dd98:	mov    rbx,rax
  94dd9b:	call   920350 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  94dda0:	mov    rdi,rbx
  94dda3:	call   e1b80 <_Unwind_Resume@plt>

Disassembly of section .fini:
