; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000949b30 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long):
  949b30:	endbr64
  949b34:	push   r15
  949b36:	push   r14
  949b38:	push   r13
  949b3a:	mov    r13,rsi
  949b3d:	push   r12
  949b3f:	mov    r12,rcx
  949b42:	push   rbp
  949b43:	mov    rbp,rdi
  949b46:	push   rbx
  949b47:	mov    rbx,rdx
  949b4a:	sub    rsp,0x148
  949b51:	mov    rax,QWORD PTR fs:0x28
  949b5a:	mov    QWORD PTR [rsp+0x138],rax
  949b62:	xor    eax,eax
  949b64:	movzx  edx,BYTE PTR [rip+0x2d5572]        # c1f0dd <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  949b6b:	test   dl,dl
  949b6d:	je     94a27a <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x74a>
  949b73:	mov    BYTE PTR [rsp],dl
  949b76:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  949b82:	cmp    QWORD PTR [rax-0x1d8],0x0
  949b8a:	mov    r14,rax
  949b8d:	lea    rcx,[rax-0x360]
  949b94:	je     949bbf <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x8f>
  949b96:	movzx  r15d,BYTE PTR [rax-0x1e0]
  949b9e:	cmp    r15b,0xf
  949ba2:	jbe    94a2dd <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x7ad>
  949ba8:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  949bb4:	add    r15d,0x1
  949bb8:	mov    BYTE PTR [rax-0x1e0],r15b
  949bbf:	mov    rdx,QWORD PTR [r13+0x38]
  949bc3:	mov    rcx,QWORD PTR [r13+0x48]
  949bc7:	mov    r9,QWORD PTR [rbp+0x8]
  949bcb:	mov    rdi,QWORD PTR [rbp+0x0]
  949bcf:	mov    r8,rdx
  949bd2:	imul   r8,rcx
  949bd6:	mov    rax,r9
  949bd9:	sub    rax,rdi
  949bdc:	sar    rax,0x2
  949be0:	lea    rsi,[r8*4+0x0]
  949be8:	cmp    rax,rsi
  949beb:	jb     94a2c5 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x795>
  949bf1:	cmp    rsi,rax
  949bf4:	jb     94a261 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x731>
  949bfa:	mov    rax,QWORD PTR [r13+0x40]
  949bfe:	lea    r9,[rcx+rcx*1]
  949c02:	imul   rdx,r9
  949c06:	mov    QWORD PTR [rsp+0xa0],rax
  949c0e:	test   rax,rax
  949c11:	je     94a103 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x5d3>
  949c17:	mov    r10,QWORD PTR [r13+0x0]
  949c1b:	mov    rdi,QWORD PTR [rbx]
  949c1e:	lea    r11,[r9+rcx*1]
  949c22:	xor    r15d,r15d
  949c25:	mov    rsi,QWORD PTR [rbx+0x10]
  949c29:	mov    rbx,rcx
  949c2c:	mov    r8,QWORD PTR [rbp+0x0]
  949c30:	shl    r11,0x2
  949c34:	shl    rbx,0x4
  949c38:	lea    r14,[rcx*4+0x0]
  949c40:	lea    rcx,[r10+rbx*1]
  949c44:	lea    rax,[r14+rdx*1]
  949c48:	mov    QWORD PTR [rsp+0x18],rbx
  949c4d:	mov    QWORD PTR [rsp+0xd8],rcx
  949c55:	lea    rcx,[r8+rbx*1]
  949c59:	mov    QWORD PTR [rsp+0xd0],rcx
  949c61:	lea    rcx,[rdx*4+0x0]
  949c69:	lea    rdx,[r8+rcx*1]
  949c6d:	mov    QWORD PTR [rsp+0xa8],r14
  949c75:	mov    QWORD PTR [rsp+0xc0],rdx
  949c7d:	mov    rdx,rax
  949c80:	lea    rax,[r8+rax*4]
  949c84:	mov    QWORD PTR [rsp+0xb0],rax
  949c8c:	lea    rax,[rsp+0x108]
  949c94:	sub    rdx,r9
  949c97:	mov    QWORD PTR [rsp+0x20],rax
  949c9c:	lea    rax,[rsp+0x110]
  949ca4:	mov    QWORD PTR [rsp+0x28],rax
  949ca9:	lea    rax,[rsp+0x118]
  949cb1:	mov    QWORD PTR [rsp+0x30],rax
  949cb6:	lea    rax,[rsp+0x120]
  949cbe:	mov    QWORD PTR [rsp+0x38],rax
  949cc3:	lea    rax,[rsp+0x128]
  949ccb:	mov    QWORD PTR [rsp+0xc8],rcx
  949cd3:	lea    rcx,[r8+rdx*4]
  949cd7:	xor    edx,edx
  949cd9:	mov    QWORD PTR [rsp+0x40],rax
  949cde:	lea    rax,[rsp+0x130]
  949ce6:	mov    QWORD PTR [rsp+0xb8],rcx
  949cee:	mov    rcx,r14
  949cf1:	mov    QWORD PTR [rsp+0x48],rax
  949cf6:	mov    rax,r12
  949cf9:	mov    r12,r15
  949cfc:	nop    DWORD PTR [rax+0x0]
  949d00:	movss  xmm0,DWORD PTR [r10+rdx*4]
  949d06:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949d0b:	movaps xmm4,xmm0
  949d0e:	movss  xmm0,DWORD PTR [r10+rcx*1]
  949d14:	shufps xmm4,xmm4,0xe0
  949d18:	movq   rbp,xmm4
  949d1d:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949d22:	movaps xmm15,xmm0
  949d26:	movss  xmm0,DWORD PTR [r10+r9*4]
  949d2c:	shufps xmm15,xmm15,0xe0
  949d31:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949d36:	movaps xmm10,xmm0
  949d3a:	movss  xmm0,DWORD PTR [r10+r11*1]
  949d40:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949d45:	mov    rbx,QWORD PTR [rsp+0xd8]
  949d4d:	movaps xmm9,xmm0
  949d51:	movss  xmm0,DWORD PTR [rbx+rdx*4]
  949d56:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949d5b:	movaps xmm14,xmm0
  949d5f:	movss  xmm0,DWORD PTR [rbx+rcx*1]
  949d64:	shufps xmm14,xmm14,0xe0
  949d69:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949d6e:	movaps xmm13,xmm0
  949d72:	movss  xmm0,DWORD PTR [r8+rdx*4]
  949d78:	shufps xmm13,xmm13,0xe0
  949d7d:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949d82:	movaps xmm7,xmm0
  949d85:	movss  xmm0,DWORD PTR [r8+rcx*1]
  949d8b:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949d90:	movaps xmm12,xmm0
  949d94:	movss  xmm0,DWORD PTR [r8+r9*4]
  949d9a:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949d9f:	movaps xmm5,xmm0
  949da2:	movss  xmm0,DWORD PTR [r8+r11*1]
  949da8:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949dad:	mov    rbx,QWORD PTR [rsp+0xd0]
  949db5:	movaps xmm6,xmm0
  949db8:	movss  xmm0,DWORD PTR [rbx+rdx*4]
  949dbd:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949dc2:	movaps xmm2,xmm0
  949dc5:	movss  xmm0,DWORD PTR [rbx+rcx*1]
  949dca:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949dcf:	mov    rbx,QWORD PTR [rsp+0xc0]
  949dd7:	movaps xmm8,xmm0
  949ddb:	movss  xmm0,DWORD PTR [rbx+rdx*4]
  949de0:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949de5:	unpcklps xmm7,xmm0
  949de8:	movss  xmm0,DWORD PTR [rbx+rcx*1]
  949ded:	movlps QWORD PTR [rsp+0x108],xmm7
  949df5:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949dfa:	unpcklps xmm12,xmm0
  949dfe:	movlps QWORD PTR [rsp+0x110],xmm12
  949e07:	mov    rbx,QWORD PTR [rsp+0xb8]
  949e0f:	movss  xmm0,DWORD PTR [rbx+rdx*4]
  949e14:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949e19:	movaps xmm3,xmm0
  949e1c:	movaps xmm0,xmm5
  949e1f:	unpcklps xmm0,xmm3
  949e22:	movlps QWORD PTR [rsp+0x118],xmm0
  949e2a:	movss  xmm0,DWORD PTR [rbx+rcx*1]
  949e2f:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949e34:	mov    rbx,QWORD PTR [rsp+0xb0]
  949e3c:	movaps xmm4,xmm0
  949e3f:	movaps xmm0,xmm6
  949e42:	unpcklps xmm0,xmm4
  949e45:	movlps QWORD PTR [rsp+0x120],xmm0
  949e4d:	movss  xmm0,DWORD PTR [rbx+rdx*4]
  949e52:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949e57:	unpcklps xmm2,xmm0
  949e5a:	movss  xmm0,DWORD PTR [rbx+rcx*1]
  949e5f:	xor    ebx,ebx
  949e61:	movlps QWORD PTR [rsp+0x128],xmm2
  949e69:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949e6e:	unpcklps xmm8,xmm0
  949e72:	movlps QWORD PTR [rsp+0x130],xmm8
  949e7b:	test   rax,rax
  949e7e:	je     949f9c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x46c>
  949e84:	nop    DWORD PTR [rax+0x0]
  949e88:	movaps xmm0,xmm8
  949e8c:	movaps xmm8,xmm2
  949e90:	movaps xmm11,xmm10
  949e94:	mulps  xmm8,xmm13
  949e98:	movaps xmm1,xmm0
  949e9b:	mulps  xmm1,xmm14
  949e9f:	mulps  xmm0,xmm13
  949ea3:	mulps  xmm2,xmm14
  949ea7:	mulss  xmm11,xmm5
  949eac:	mulss  xmm5,xmm9
  949eb1:	addps  xmm8,xmm1
  949eb5:	movq   xmm1,rbp
  949eba:	subps  xmm2,xmm0
  949ebd:	movaps xmm0,xmm12
  949ec1:	movaps xmm12,xmm7
  949ec5:	mulps  xmm1,xmm0
  949ec8:	mulps  xmm12,xmm15
  949ecc:	mulps  xmm0,xmm15
  949ed0:	addps  xmm12,xmm1
  949ed4:	movq   xmm1,rbp
  949ed9:	mulps  xmm7,xmm1
  949edc:	subps  xmm7,xmm0
  949edf:	movss  xmm0,DWORD PTR [rdi+rbx*4]
  949ee4:	call   9205f0 <float cpl::simd::broadcast<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949ee9:	movaps xmm1,xmm0
  949eec:	movaps xmm0,xmm9
  949ef0:	mulss  xmm0,xmm6
  949ef4:	mulss  xmm6,xmm10
  949ef9:	subss  xmm11,xmm0
  949efe:	movss  xmm0,DWORD PTR [rsi+rbx*4]
  949f03:	add    rbx,0x1
  949f07:	addss  xmm6,xmm5
  949f0b:	call   9205f0 <float cpl::simd::broadcast<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949f10:	movaps xmm5,xmm11
  949f14:	movaps xmm11,xmm10
  949f18:	addss  xmm5,xmm1
  949f1c:	mulss  xmm11,xmm3
  949f21:	unpcklps xmm1,xmm0
  949f24:	addps  xmm2,xmm1
  949f27:	addps  xmm7,xmm1
  949f2a:	movaps xmm1,xmm4
  949f2d:	mulss  xmm1,xmm9
  949f32:	mulss  xmm3,xmm9
  949f37:	mulss  xmm4,xmm10
  949f3c:	subss  xmm11,xmm1
  949f41:	addss  xmm4,xmm3
  949f45:	movaps xmm3,xmm0
  949f48:	addss  xmm3,xmm11
  949f4d:	cmp    rax,rbx
  949f50:	jne    949e88 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x358>
  949f56:	movss  DWORD PTR [rsp+0x120],xmm6
  949f5f:	movss  DWORD PTR [rsp+0x118],xmm5
  949f68:	movlps QWORD PTR [rsp+0x110],xmm12
  949f71:	movlps QWORD PTR [rsp+0x108],xmm7
  949f79:	movss  DWORD PTR [rsp+0x124],xmm4
  949f82:	movss  DWORD PTR [rsp+0x11c],xmm3
  949f8b:	movlps QWORD PTR [rsp+0x130],xmm8
  949f94:	movlps QWORD PTR [rsp+0x128],xmm2
  949f9c:	mov    rbx,r9
  949f9f:	mov    r13,r8
  949fa2:	mov    QWORD PTR [rsp+0x50],rdx
  949fa7:	lea    rbp,[r8+rcx*1]
  949fab:	sub    rbx,rdx
  949fae:	mov    QWORD PTR [rsp+0x60],r9
  949fb3:	sub    r13,r12
  949fb6:	shl    rbx,0x2
  949fba:	mov    QWORD PTR [rsp+0x68],r10
  949fbf:	mov    QWORD PTR [rsp],rbx
  949fc3:	mov    rbx,QWORD PTR [rsp+0xa8]
  949fcb:	mov    QWORD PTR [rsp+0x70],rdi
  949fd0:	lea    r15,[rbx+rcx*1]
  949fd4:	lea    r14,[rbx+r11*1]
  949fd8:	mov    QWORD PTR [rsp+0x78],r8
  949fdd:	xor    ebx,ebx
  949fdf:	mov    QWORD PTR [rsp+0x80],rsi
  949fe7:	add    r15,r12
  949fea:	add    r14,r12
  949fed:	mov    QWORD PTR [rsp+0x88],r11
  949ff5:	mov    QWORD PTR [rsp+0x90],rcx
  949ffd:	mov    QWORD PTR [rsp+0x98],rax
  94a005:	mov    QWORD PTR [rsp+0x58],r12
  94a00a:	mov    r12,QWORD PTR [rsp+0xc8]
  94a012:	mov    rax,QWORD PTR [rsp+0x20]
  94a017:	mov    rdi,r13
  94a01a:	movss  xmm0,DWORD PTR [rax+rbx*4]
  94a01f:	call   85fb30 <cpl::simd::store(float*, float)>
  94a024:	mov    rax,QWORD PTR [rsp+0x28]
  94a029:	mov    rdi,rbp
  94a02c:	movss  xmm0,DWORD PTR [rax+rbx*4]
  94a031:	call   85fb30 <cpl::simd::store(float*, float)>
  94a036:	mov    rax,QWORD PTR [rsp+0x30]
  94a03b:	movss  xmm0,DWORD PTR [rax+rbx*4]
  94a040:	mov    rax,QWORD PTR [rsp]
  94a044:	lea    rdi,[rax+r13*1]
  94a048:	call   85fb30 <cpl::simd::store(float*, float)>
  94a04d:	mov    rax,QWORD PTR [rsp+0x38]
  94a052:	lea    rdi,[r15+rbp*1]
  94a056:	movss  xmm0,DWORD PTR [rax+rbx*4]
  94a05b:	call   85fb30 <cpl::simd::store(float*, float)>
  94a060:	mov    rax,QWORD PTR [rsp+0x40]
  94a065:	movss  xmm0,DWORD PTR [rax+rbx*4]
  94a06a:	mov    rax,QWORD PTR [rsp+0x18]
  94a06f:	lea    rdi,[rax+r13*1]
  94a073:	add    r13,r12
  94a076:	call   85fb30 <cpl::simd::store(float*, float)>
  94a07b:	mov    rax,QWORD PTR [rsp+0x48]
  94a080:	lea    rdi,[r14+rbp*1]
  94a084:	add    rbp,r12
  94a087:	movss  xmm0,DWORD PTR [rax+rbx*4]
  94a08c:	call   85fb30 <cpl::simd::store(float*, float)>
  94a091:	test   rbx,rbx
  94a094:	jne    94a0a0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x570>
  94a096:	mov    ebx,0x1
  94a09b:	jmp    94a012 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4e2>
  94a0a0:	mov    rdx,QWORD PTR [rsp+0x50]
  94a0a5:	mov    r12,QWORD PTR [rsp+0x58]
  94a0aa:	mov    r9,QWORD PTR [rsp+0x60]
  94a0af:	mov    r11,QWORD PTR [rsp+0x88]
  94a0b7:	mov    rcx,QWORD PTR [rsp+0x90]
  94a0bf:	add    rdx,0x1
  94a0c3:	mov    r10,QWORD PTR [rsp+0x68]
  94a0c8:	sub    r12,0x4
  94a0cc:	mov    rdi,QWORD PTR [rsp+0x70]
  94a0d1:	mov    rbx,QWORD PTR [rsp+0xa0]
  94a0d9:	add    r9,0x1
  94a0dd:	add    r11,0x4
  94a0e1:	mov    r8,QWORD PTR [rsp+0x78]
  94a0e6:	mov    rsi,QWORD PTR [rsp+0x80]
  94a0ee:	add    rcx,0x4
  94a0f2:	mov    rax,QWORD PTR [rsp+0x98]
  94a0fa:	cmp    rdx,rbx
  94a0fd:	jne    949d00 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x1d0>
  94a103:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94a10f:	lea    rbp,[rax-0x360]
  94a116:	mov    r12,rax
  94a119:	mov    r13,QWORD PTR [rbp+0x188]
  94a120:	test   r13,r13
  94a123:	je     94a143 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x613>
  94a125:	movzx  eax,BYTE PTR [rbp+0x180]
  94a12c:	lea    ebx,[rax-0x1]
  94a12f:	mov    BYTE PTR [rbp+0x180],bl
  94a135:	cmp    bl,0xf
  94a138:	jbe    94a16c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x63c>
  94a13a:	add    WORD PTR [r13+0xfe2],0x1
  94a143:	mov    rax,QWORD PTR [rsp+0x138]
  94a14b:	sub    rax,QWORD PTR fs:0x28
  94a154:	jne    94a370 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x840>
  94a15a:	add    rsp,0x148
  94a161:	pop    rbx
  94a162:	pop    rbp
  94a163:	pop    r12
  94a165:	pop    r13
  94a167:	pop    r14
  94a169:	pop    r15
  94a16b:	ret
  94a16c:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  94a171:	movzx  r15d,bl
  94a175:	lea    rdx,[r15+r15*2]
  94a179:	mov    rcx,rax
  94a17c:	movq   xmm0,rax
  94a181:	shl    rdx,0x3
  94a185:	movdqu xmm1,XMMWORD PTR [r12+rdx*1-0x360]
  94a18f:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  94a194:	movq   xmm6,rcx
  94a199:	movzx  r12d,BYTE PTR [rbp+0x181]
  94a1a1:	punpcklqdq xmm0,xmm6
  94a1a5:	psubq  xmm0,xmm1
  94a1a9:	cmp    r12b,bl
  94a1ac:	jae    94a1cc <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x69c>
  94a1ae:	lea    eax,[r15-0x1]
  94a1b2:	movdqa xmm1,xmm0
  94a1b6:	cdqe
  94a1b8:	lea    rax,[rax+rax*2]
  94a1bc:	movq   xmm2,QWORD PTR [rbp+rax*8+0x8]
  94a1c2:	paddq  xmm1,xmm2
  94a1c6:	movq   QWORD PTR [rbp+rax*8+0x8],xmm1
  94a1cc:	movzx  ebp,WORD PTR [r13+0xfe0]
  94a1d4:	cmp    bp,0x7f
  94a1d8:	je     94a336 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x806>
  94a1de:	movaps XMMWORD PTR [rsp],xmm0
  94a1e2:	sub    ebx,r12d
  94a1e5:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94a1f1:	lea    rdx,[r15+r15*2]
  94a1f5:	movdqa xmm0,XMMWORD PTR [rsp]
  94a1fa:	mov    DWORD PTR [rsp+0xf8],0x0
  94a205:	movups XMMWORD PTR [rsp+0xe8],xmm0
  94a20d:	lea    rdx,[rax+rdx*8-0x360]
  94a215:	mov    rax,QWORD PTR [rdx]
  94a218:	mov    QWORD PTR [rsp+0xe0],rax
  94a220:	movzx  eax,BYTE PTR [rdx+0x10]
  94a224:	lea    edx,[rbp+0x1]
  94a227:	movdqa xmm6,XMMWORD PTR [rsp+0xe0]
  94a230:	mov    WORD PTR [r13+0xfe0],dx
  94a238:	mov    ah,bl
  94a23a:	mov    WORD PTR [rsp+0xfc],ax
  94a242:	movzx  eax,bp
  94a245:	shl    rax,0x5
  94a249:	add    rax,r13
  94a24c:	movups XMMWORD PTR [rax],xmm6
  94a24f:	movdqu xmm6,XMMWORD PTR [rsp+0xee]
  94a258:	movups XMMWORD PTR [rax+0xe],xmm6
  94a25c:	jmp    94a143 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x613>
  94a261:	shl    r8,0x4
  94a265:	add    rdi,r8
  94a268:	cmp    r9,rdi
  94a26b:	je     949bfa <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xca>
  94a271:	mov    QWORD PTR [rbp+0x8],rdi
  94a275:	jmp    949bfa <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xca>
  94a27a:	mov    eax,0x1
  94a27f:	lock xadd DWORD PTR [rip+0x2c98d9],eax        # c13b60 <cpl::Profiling::registerRegion(char const*)::counter>
  94a287:	mov    ecx,0x1
  94a28c:	mov    edx,0x1
  94a291:	add    eax,0x2
  94a294:	cmp    eax,0xfe
  94a299:	ja     94a2b3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x783>
  94a29b:	mov    ecx,eax
  94a29d:	lea    rdx,[rip+0x2c7d7c]        # c12020 <cpl::Profiling::regions>
  94a2a4:	lea    rsi,[rip+0x9de9b]        # 9e8146 <_IO_stdin_used+0x18146>
  94a2ab:	mov    QWORD PTR [rdx+rcx*8],rsi
  94a2af:	mov    edx,eax
  94a2b1:	mov    ecx,eax
  94a2b3:	xor    eax,eax
  94a2b5:	lock cmpxchg BYTE PTR [rip+0x2d4e20],cl        # c1f0dd <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  94a2bd:	cmovne edx,eax
  94a2c0:	jmp    949b73 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x43>
  94a2c5:	sub    rsi,rax
  94a2c8:	mov    rdi,rbp
  94a2cb:	call   91d820 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  94a2d0:	mov    rdx,QWORD PTR [r13+0x38]
  94a2d4:	mov    rcx,QWORD PTR [r13+0x48]
  94a2d8:	jmp    949bfa <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xca>
  94a2dd:	mov    QWORD PTR [rsp+0x18],rcx
  94a2e2:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  94a2e7:	movzx  edx,BYTE PTR [rsp]
  94a2eb:	mov    rcx,QWORD PTR [rsp+0x18]
  94a2f0:	mov    QWORD PTR [rsp+0xe8],0x0
  94a2fc:	mov    QWORD PTR [rsp+0xe0],rax
  94a304:	lea    rax,[r15+r15*2]
  94a308:	movdqa xmm6,XMMWORD PTR [rsp+0xe0]
  94a311:	shl    rax,0x3
  94a315:	mov    BYTE PTR [rsp+0xf0],dl
  94a31c:	movups XMMWORD PTR [r14+rax*1-0x360],xmm6
  94a325:	mov    BYTE PTR [rax+rcx*1+0x10],dl
  94a329:	movzx  r15d,BYTE PTR [rcx+0x180]
  94a331:	jmp    949ba8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x78>
  94a336:	movzx  eax,WORD PTR [r13+0xfe2]
  94a33e:	cmp    ax,0xffff
  94a342:	je     94a143 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x613>
  94a348:	add    eax,0x1
  94a34b:	mov    WORD PTR [r13+0xfe2],ax
  94a353:	jmp    94a143 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x613>
  94a358:	call   920350 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  94a35d:	mov    rax,QWORD PTR [rsp+0x138]
  94a365:	sub    rax,QWORD PTR fs:0x28
  94a36e:	je     94a37e <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x84e>
  94a370:	call   e0f90 <__stack_chk_fail@plt>
  94a375:	endbr64
  94a379:	mov    rbx,rax
  94a37c:	jmp    94a358 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x828>
  94a37e:	mov    rdi,rbx
  94a381:	call   e1b80 <_Unwind_Resume@plt>

Disassembly of section .fini:
