; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

000000000094ddb0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long):
  94ddb0:	endbr64
  94ddb4:	push   r15
  94ddb6:	push   r14
  94ddb8:	push   r13
  94ddba:	mov    r13,rsi
  94ddbd:	push   r12
  94ddbf:	mov    r12,rcx
  94ddc2:	push   rbp
  94ddc3:	mov    rbp,rdi
  94ddc6:	push   rbx
  94ddc7:	mov    rbx,rdx
  94ddca:	sub    rsp,0x148
  94ddd1:	mov    rax,QWORD PTR fs:0x28
  94ddda:	mov    QWORD PTR [rsp+0x138],rax
  94dde2:	xor    eax,eax
  94dde4:	movzx  edx,BYTE PTR [rip+0x2d12e8]        # c1f0d3 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  94ddeb:	test   dl,dl
  94dded:	je     94e4fa <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x74a>
  94ddf3:	mov    BYTE PTR [rsp],dl
  94ddf6:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94de02:	cmp    QWORD PTR [rax-0x1d8],0x0
  94de0a:	mov    r14,rax
  94de0d:	lea    rcx,[rax-0x360]
  94de14:	je     94de3f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x8f>
  94de16:	movzx  r15d,BYTE PTR [rax-0x1e0]
  94de1e:	cmp    r15b,0xf
  94de22:	jbe    94e55d <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x7ad>
  94de28:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94de34:	add    r15d,0x1
  94de38:	mov    BYTE PTR [rax-0x1e0],r15b
  94de3f:	mov    rdx,QWORD PTR [r13+0x38]
  94de43:	mov    rcx,QWORD PTR [r13+0x48]
  94de47:	mov    r9,QWORD PTR [rbp+0x8]
  94de4b:	mov    rdi,QWORD PTR [rbp+0x0]
  94de4f:	mov    r8,rdx
  94de52:	imul   r8,rcx
  94de56:	mov    rax,r9
  94de59:	sub    rax,rdi
  94de5c:	sar    rax,0x2
  94de60:	lea    rsi,[r8*4+0x0]
  94de68:	cmp    rax,rsi
  94de6b:	jb     94e545 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x795>
  94de71:	cmp    rsi,rax
  94de74:	jb     94e4e1 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x731>
  94de7a:	mov    rax,QWORD PTR [r13+0x40]
  94de7e:	lea    r9,[rcx+rcx*1]
  94de82:	imul   rdx,r9
  94de86:	mov    QWORD PTR [rsp+0xa0],rax
  94de8e:	test   rax,rax
  94de91:	je     94e383 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x5d3>
  94de97:	mov    r10,QWORD PTR [r13+0x0]
  94de9b:	mov    rdi,QWORD PTR [rbx]
  94de9e:	lea    r11,[r9+rcx*1]
  94dea2:	xor    r15d,r15d
  94dea5:	mov    rsi,QWORD PTR [rbx+0x10]
  94dea9:	mov    rbx,rcx
  94deac:	mov    r8,QWORD PTR [rbp+0x0]
  94deb0:	shl    r11,0x2
  94deb4:	shl    rbx,0x4
  94deb8:	lea    r14,[rcx*4+0x0]
  94dec0:	lea    rcx,[r10+rbx*1]
  94dec4:	lea    rax,[r14+rdx*1]
  94dec8:	mov    QWORD PTR [rsp+0x18],rbx
  94decd:	mov    QWORD PTR [rsp+0xd8],rcx
  94ded5:	lea    rcx,[r8+rbx*1]
  94ded9:	mov    QWORD PTR [rsp+0xd0],rcx
  94dee1:	lea    rcx,[rdx*4+0x0]
  94dee9:	lea    rdx,[r8+rcx*1]
  94deed:	mov    QWORD PTR [rsp+0xa8],r14
  94def5:	mov    QWORD PTR [rsp+0xc0],rdx
  94defd:	mov    rdx,rax
  94df00:	lea    rax,[r8+rax*4]
  94df04:	mov    QWORD PTR [rsp+0xb0],rax
  94df0c:	lea    rax,[rsp+0x108]
  94df14:	sub    rdx,r9
  94df17:	mov    QWORD PTR [rsp+0x20],rax
  94df1c:	lea    rax,[rsp+0x110]
  94df24:	mov    QWORD PTR [rsp+0x28],rax
  94df29:	lea    rax,[rsp+0x118]
  94df31:	mov    QWORD PTR [rsp+0x30],rax
  94df36:	lea    rax,[rsp+0x120]
  94df3e:	mov    QWORD PTR [rsp+0x38],rax
  94df43:	lea    rax,[rsp+0x128]
  94df4b:	mov    QWORD PTR [rsp+0xc8],rcx
  94df53:	lea    rcx,[r8+rdx*4]
  94df57:	xor    edx,edx
  94df59:	mov    QWORD PTR [rsp+0x40],rax
  94df5e:	lea    rax,[rsp+0x130]
  94df66:	mov    QWORD PTR [rsp+0xb8],rcx
  94df6e:	mov    rcx,r14
  94df71:	mov    QWORD PTR [rsp+0x48],rax
  94df76:	mov    rax,r12
  94df79:	mov    r12,r15
  94df7c:	nop    DWORD PTR [rax+0x0]
  94df80:	movss  xmm0,DWORD PTR [r10+rdx*4]
  94df86:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94df8b:	movaps xmm4,xmm0
  94df8e:	movss  xmm0,DWORD PTR [r10+rcx*1]
  94df94:	shufps xmm4,xmm4,0xe0
  94df98:	movq   rbp,xmm4
  94df9d:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94dfa2:	movaps xmm15,xmm0
  94dfa6:	movss  xmm0,DWORD PTR [r10+r9*4]
  94dfac:	shufps xmm15,xmm15,0xe0
  94dfb1:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94dfb6:	movaps xmm10,xmm0
  94dfba:	movss  xmm0,DWORD PTR [r10+r11*1]
  94dfc0:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94dfc5:	mov    rbx,QWORD PTR [rsp+0xd8]
  94dfcd:	movaps xmm9,xmm0
  94dfd1:	movss  xmm0,DWORD PTR [rbx+rdx*4]
  94dfd6:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94dfdb:	movaps xmm14,xmm0
  94dfdf:	movss  xmm0,DWORD PTR [rbx+rcx*1]
  94dfe4:	shufps xmm14,xmm14,0xe0
  94dfe9:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94dfee:	movaps xmm13,xmm0
  94dff2:	movss  xmm0,DWORD PTR [r8+rdx*4]
  94dff8:	shufps xmm13,xmm13,0xe0
  94dffd:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e002:	movaps xmm7,xmm0
  94e005:	movss  xmm0,DWORD PTR [r8+rcx*1]
  94e00b:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e010:	movaps xmm12,xmm0
  94e014:	movss  xmm0,DWORD PTR [r8+r9*4]
  94e01a:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e01f:	movaps xmm5,xmm0
  94e022:	movss  xmm0,DWORD PTR [r8+r11*1]
  94e028:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e02d:	mov    rbx,QWORD PTR [rsp+0xd0]
  94e035:	movaps xmm6,xmm0
  94e038:	movss  xmm0,DWORD PTR [rbx+rdx*4]
  94e03d:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e042:	movaps xmm2,xmm0
  94e045:	movss  xmm0,DWORD PTR [rbx+rcx*1]
  94e04a:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e04f:	mov    rbx,QWORD PTR [rsp+0xc0]
  94e057:	movaps xmm8,xmm0
  94e05b:	movss  xmm0,DWORD PTR [rbx+rdx*4]
  94e060:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e065:	unpcklps xmm7,xmm0
  94e068:	movss  xmm0,DWORD PTR [rbx+rcx*1]
  94e06d:	movlps QWORD PTR [rsp+0x108],xmm7
  94e075:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e07a:	unpcklps xmm12,xmm0
  94e07e:	movlps QWORD PTR [rsp+0x110],xmm12
  94e087:	mov    rbx,QWORD PTR [rsp+0xb8]
  94e08f:	movss  xmm0,DWORD PTR [rbx+rdx*4]
  94e094:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e099:	movaps xmm3,xmm0
  94e09c:	movaps xmm0,xmm5
  94e09f:	unpcklps xmm0,xmm3
  94e0a2:	movlps QWORD PTR [rsp+0x118],xmm0
  94e0aa:	movss  xmm0,DWORD PTR [rbx+rcx*1]
  94e0af:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e0b4:	mov    rbx,QWORD PTR [rsp+0xb0]
  94e0bc:	movaps xmm4,xmm0
  94e0bf:	movaps xmm0,xmm6
  94e0c2:	unpcklps xmm0,xmm4
  94e0c5:	movlps QWORD PTR [rsp+0x120],xmm0
  94e0cd:	movss  xmm0,DWORD PTR [rbx+rdx*4]
  94e0d2:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e0d7:	unpcklps xmm2,xmm0
  94e0da:	movss  xmm0,DWORD PTR [rbx+rcx*1]
  94e0df:	xor    ebx,ebx
  94e0e1:	movlps QWORD PTR [rsp+0x128],xmm2
  94e0e9:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e0ee:	unpcklps xmm8,xmm0
  94e0f2:	movlps QWORD PTR [rsp+0x130],xmm8
  94e0fb:	test   rax,rax
  94e0fe:	je     94e21c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x46c>
  94e104:	nop    DWORD PTR [rax+0x0]
  94e108:	movaps xmm0,xmm8
  94e10c:	movaps xmm8,xmm2
  94e110:	movaps xmm11,xmm10
  94e114:	mulps  xmm8,xmm13
  94e118:	movaps xmm1,xmm0
  94e11b:	mulps  xmm1,xmm14
  94e11f:	mulps  xmm0,xmm13
  94e123:	mulps  xmm2,xmm14
  94e127:	mulss  xmm11,xmm5
  94e12c:	mulss  xmm5,xmm9
  94e131:	addps  xmm8,xmm1
  94e135:	movq   xmm1,rbp
  94e13a:	subps  xmm2,xmm0
  94e13d:	movaps xmm0,xmm12
  94e141:	movaps xmm12,xmm7
  94e145:	mulps  xmm1,xmm0
  94e148:	mulps  xmm12,xmm15
  94e14c:	mulps  xmm0,xmm15
  94e150:	addps  xmm12,xmm1
  94e154:	movq   xmm1,rbp
  94e159:	mulps  xmm7,xmm1
  94e15c:	subps  xmm7,xmm0
  94e15f:	movss  xmm0,DWORD PTR [rdi+rbx*4]
  94e164:	call   9205f0 <float cpl::simd::broadcast<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e169:	movaps xmm1,xmm0
  94e16c:	movaps xmm0,xmm9
  94e170:	mulss  xmm0,xmm6
  94e174:	mulss  xmm6,xmm10
  94e179:	subss  xmm11,xmm0
  94e17e:	movss  xmm0,DWORD PTR [rsi+rbx*4]
  94e183:	add    rbx,0x1
  94e187:	addss  xmm6,xmm5
  94e18b:	call   9205f0 <float cpl::simd::broadcast<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94e190:	movaps xmm5,xmm11
  94e194:	movaps xmm11,xmm10
  94e198:	addss  xmm5,xmm1
  94e19c:	mulss  xmm11,xmm3
  94e1a1:	unpcklps xmm1,xmm0
  94e1a4:	addps  xmm2,xmm1
  94e1a7:	addps  xmm7,xmm1
  94e1aa:	movaps xmm1,xmm4
  94e1ad:	mulss  xmm1,xmm9
  94e1b2:	mulss  xmm3,xmm9
  94e1b7:	mulss  xmm4,xmm10
  94e1bc:	subss  xmm11,xmm1
  94e1c1:	addss  xmm4,xmm3
  94e1c5:	movaps xmm3,xmm0
  94e1c8:	addss  xmm3,xmm11
  94e1cd:	cmp    rax,rbx
  94e1d0:	jne    94e108 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x358>
  94e1d6:	movss  DWORD PTR [rsp+0x120],xmm6
  94e1df:	movss  DWORD PTR [rsp+0x118],xmm5
  94e1e8:	movlps QWORD PTR [rsp+0x110],xmm12
  94e1f1:	movlps QWORD PTR [rsp+0x108],xmm7
  94e1f9:	movss  DWORD PTR [rsp+0x124],xmm4
  94e202:	movss  DWORD PTR [rsp+0x11c],xmm3
  94e20b:	movlps QWORD PTR [rsp+0x130],xmm8
  94e214:	movlps QWORD PTR [rsp+0x128],xmm2
  94e21c:	mov    rbx,r9
  94e21f:	mov    r13,r8
  94e222:	mov    QWORD PTR [rsp+0x50],rdx
  94e227:	lea    rbp,[r8+rcx*1]
  94e22b:	sub    rbx,rdx
  94e22e:	mov    QWORD PTR [rsp+0x60],r9
  94e233:	sub    r13,r12
  94e236:	shl    rbx,0x2
  94e23a:	mov    QWORD PTR [rsp+0x68],r10
  94e23f:	mov    QWORD PTR [rsp],rbx
  94e243:	mov    rbx,QWORD PTR [rsp+0xa8]
  94e24b:	mov    QWORD PTR [rsp+0x70],rdi
  94e250:	lea    r15,[rbx+rcx*1]
  94e254:	lea    r14,[rbx+r11*1]
  94e258:	mov    QWORD PTR [rsp+0x78],r8
  94e25d:	xor    ebx,ebx
  94e25f:	mov    QWORD PTR [rsp+0x80],rsi
  94e267:	add    r15,r12
  94e26a:	add    r14,r12
  94e26d:	mov    QWORD PTR [rsp+0x88],r11
  94e275:	mov    QWORD PTR [rsp+0x90],rcx
  94e27d:	mov    QWORD PTR [rsp+0x98],rax
  94e285:	mov    QWORD PTR [rsp+0x58],r12
  94e28a:	mov    r12,QWORD PTR [rsp+0xc8]
  94e292:	mov    rax,QWORD PTR [rsp+0x20]
  94e297:	mov    rdi,r13
  94e29a:	movss  xmm0,DWORD PTR [rax+rbx*4]
  94e29f:	call   85fb30 <cpl::simd::store(float*, float)>
  94e2a4:	mov    rax,QWORD PTR [rsp+0x28]
  94e2a9:	mov    rdi,rbp
  94e2ac:	movss  xmm0,DWORD PTR [rax+rbx*4]
  94e2b1:	call   85fb30 <cpl::simd::store(float*, float)>
  94e2b6:	mov    rax,QWORD PTR [rsp+0x30]
  94e2bb:	movss  xmm0,DWORD PTR [rax+rbx*4]
  94e2c0:	mov    rax,QWORD PTR [rsp]
  94e2c4:	lea    rdi,[rax+r13*1]
  94e2c8:	call   85fb30 <cpl::simd::store(float*, float)>
  94e2cd:	mov    rax,QWORD PTR [rsp+0x38]
  94e2d2:	lea    rdi,[r15+rbp*1]
  94e2d6:	movss  xmm0,DWORD PTR [rax+rbx*4]
  94e2db:	call   85fb30 <cpl::simd::store(float*, float)>
  94e2e0:	mov    rax,QWORD PTR [rsp+0x40]
  94e2e5:	movss  xmm0,DWORD PTR [rax+rbx*4]
  94e2ea:	mov    rax,QWORD PTR [rsp+0x18]
  94e2ef:	lea    rdi,[rax+r13*1]
  94e2f3:	add    r13,r12
  94e2f6:	call   85fb30 <cpl::simd::store(float*, float)>
  94e2fb:	mov    rax,QWORD PTR [rsp+0x48]
  94e300:	lea    rdi,[r14+rbp*1]
  94e304:	add    rbp,r12
  94e307:	movss  xmm0,DWORD PTR [rax+rbx*4]
  94e30c:	call   85fb30 <cpl::simd::store(float*, float)>
  94e311:	test   rbx,rbx
  94e314:	jne    94e320 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x570>
  94e316:	mov    ebx,0x1
  94e31b:	jmp    94e292 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4e2>
  94e320:	mov    rdx,QWORD PTR [rsp+0x50]
  94e325:	mov    r12,QWORD PTR [rsp+0x58]
  94e32a:	mov    r9,QWORD PTR [rsp+0x60]
  94e32f:	mov    r11,QWORD PTR [rsp+0x88]
  94e337:	mov    rcx,QWORD PTR [rsp+0x90]
  94e33f:	add    rdx,0x1
  94e343:	mov    r10,QWORD PTR [rsp+0x68]
  94e348:	sub    r12,0x4
  94e34c:	mov    rdi,QWORD PTR [rsp+0x70]
  94e351:	mov    rbx,QWORD PTR [rsp+0xa0]
  94e359:	add    r9,0x1
  94e35d:	add    r11,0x4
  94e361:	mov    r8,QWORD PTR [rsp+0x78]
  94e366:	mov    rsi,QWORD PTR [rsp+0x80]
  94e36e:	add    rcx,0x4
  94e372:	mov    rax,QWORD PTR [rsp+0x98]
  94e37a:	cmp    rdx,rbx
  94e37d:	jne    94df80 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x1d0>
  94e383:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94e38f:	lea    rbp,[rax-0x360]
  94e396:	mov    r12,rax
  94e399:	mov    r13,QWORD PTR [rbp+0x188]
  94e3a0:	test   r13,r13
  94e3a3:	je     94e3c3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x613>
  94e3a5:	movzx  eax,BYTE PTR [rbp+0x180]
  94e3ac:	lea    ebx,[rax-0x1]
  94e3af:	mov    BYTE PTR [rbp+0x180],bl
  94e3b5:	cmp    bl,0xf
  94e3b8:	jbe    94e3ec <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x63c>
  94e3ba:	add    WORD PTR [r13+0xfe2],0x1
  94e3c3:	mov    rax,QWORD PTR [rsp+0x138]
  94e3cb:	sub    rax,QWORD PTR fs:0x28
  94e3d4:	jne    94e5f0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x840>
  94e3da:	add    rsp,0x148
  94e3e1:	pop    rbx
  94e3e2:	pop    rbp
  94e3e3:	pop    r12
  94e3e5:	pop    r13
  94e3e7:	pop    r14
  94e3e9:	pop    r15
  94e3eb:	ret
  94e3ec:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  94e3f1:	movzx  r15d,bl
  94e3f5:	lea    rdx,[r15+r15*2]
  94e3f9:	mov    rcx,rax
  94e3fc:	movq   xmm0,rax
  94e401:	shl    rdx,0x3
  94e405:	movdqu xmm1,XMMWORD PTR [r12+rdx*1-0x360]
  94e40f:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  94e414:	movq   xmm6,rcx
  94e419:	movzx  r12d,BYTE PTR [rbp+0x181]
  94e421:	punpcklqdq xmm0,xmm6
  94e425:	psubq  xmm0,xmm1
  94e429:	cmp    r12b,bl
  94e42c:	jae    94e44c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x69c>
  94e42e:	lea    eax,[r15-0x1]
  94e432:	movdqa xmm1,xmm0
  94e436:	cdqe
  94e438:	lea    rax,[rax+rax*2]
  94e43c:	movq   xmm2,QWORD PTR [rbp+rax*8+0x8]
  94e442:	paddq  xmm1,xmm2
  94e446:	movq   QWORD PTR [rbp+rax*8+0x8],xmm1
  94e44c:	movzx  ebp,WORD PTR [r13+0xfe0]
  94e454:	cmp    bp,0x7f
  94e458:	je     94e5b6 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x806>
  94e45e:	movaps XMMWORD PTR [rsp],xmm0
  94e462:	sub    ebx,r12d
  94e465:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94e471:	lea    rdx,[r15+r15*2]
  94e475:	movdqa xmm0,XMMWORD PTR [rsp]
  94e47a:	mov    DWORD PTR [rsp+0xf8],0x0
  94e485:	movups XMMWORD PTR [rsp+0xe8],xmm0
  94e48d:	lea    rdx,[rax+rdx*8-0x360]
  94e495:	mov    rax,QWORD PTR [rdx]
  94e498:	mov    QWORD PTR [rsp+0xe0],rax
  94e4a0:	movzx  eax,BYTE PTR [rdx+0x10]
  94e4a4:	lea    edx,[rbp+0x1]
  94e4a7:	movdqa xmm6,XMMWORD PTR [rsp+0xe0]
  94e4b0:	mov    WORD PTR [r13+0xfe0],dx
  94e4b8:	mov    ah,bl
  94e4ba:	mov    WORD PTR [rsp+0xfc],ax
  94e4c2:	movzx  eax,bp
  94e4c5:	shl    rax,0x5
  94e4c9:	add    rax,r13
  94e4cc:	movups XMMWORD PTR [rax],xmm6
  94e4cf:	movdqu xmm6,XMMWORD PTR [rsp+0xee]
  94e4d8:	movups XMMWORD PTR [rax+0xe],xmm6
  94e4dc:	jmp    94e3c3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x613>
  94e4e1:	shl    r8,0x4
  94e4e5:	add    rdi,r8
  94e4e8:	cmp    r9,rdi
  94e4eb:	je     94de7a <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xca>
  94e4f1:	mov    QWORD PTR [rbp+0x8],rdi
  94e4f5:	jmp    94de7a <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xca>
  94e4fa:	mov    eax,0x1
  94e4ff:	lock xadd DWORD PTR [rip+0x2c5659],eax        # c13b60 <cpl::Profiling::registerRegion(char const*)::counter>
  94e507:	mov    ecx,0x1
  94e50c:	mov    edx,0x1
  94e511:	add    eax,0x2
  94e514:	cmp    eax,0xfe
  94e519:	ja     94e533 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x783>
  94e51b:	mov    ecx,eax
  94e51d:	lea    rdx,[rip+0x2c3afc]        # c12020 <cpl::Profiling::regions>
  94e524:	lea    rsi,[rip+0x99c1b]        # 9e8146 <_IO_stdin_used+0x18146>
  94e52b:	mov    QWORD PTR [rdx+rcx*8],rsi
  94e52f:	mov    edx,eax
  94e531:	mov    ecx,eax
  94e533:	xor    eax,eax
  94e535:	lock cmpxchg BYTE PTR [rip+0x2d0b96],cl        # c1f0d3 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  94e53d:	cmovne edx,eax
  94e540:	jmp    94ddf3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x43>
  94e545:	sub    rsi,rax
  94e548:	mov    rdi,rbp
  94e54b:	call   91d820 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  94e550:	mov    rdx,QWORD PTR [r13+0x38]
  94e554:	mov    rcx,QWORD PTR [r13+0x48]
  94e558:	jmp    94de7a <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xca>
  94e55d:	mov    QWORD PTR [rsp+0x18],rcx
  94e562:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  94e567:	movzx  edx,BYTE PTR [rsp]
  94e56b:	mov    rcx,QWORD PTR [rsp+0x18]
  94e570:	mov    QWORD PTR [rsp+0xe8],0x0
  94e57c:	mov    QWORD PTR [rsp+0xe0],rax
  94e584:	lea    rax,[r15+r15*2]
  94e588:	movdqa xmm6,XMMWORD PTR [rsp+0xe0]
  94e591:	shl    rax,0x3
  94e595:	mov    BYTE PTR [rsp+0xf0],dl
  94e59c:	movups XMMWORD PTR [r14+rax*1-0x360],xmm6
  94e5a5:	mov    BYTE PTR [rax+rcx*1+0x10],dl
  94e5a9:	movzx  r15d,BYTE PTR [rcx+0x180]
  94e5b1:	jmp    94de28 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x78>
  94e5b6:	movzx  eax,WORD PTR [r13+0xfe2]
  94e5be:	cmp    ax,0xffff
  94e5c2:	je     94e3c3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x613>
  94e5c8:	add    eax,0x1
  94e5cb:	mov    WORD PTR [r13+0xfe2],ax
  94e5d3:	jmp    94e3c3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x613>
  94e5d8:	call   920350 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  94e5dd:	mov    rax,QWORD PTR [rsp+0x138]
  94e5e5:	sub    rax,QWORD PTR fs:0x28
  94e5ee:	je     94e5fe <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x84e>
  94e5f0:	call   e0f90 <__stack_chk_fail@plt>
  94e5f5:	endbr64
  94e5f9:	mov    rbx,rax
  94e5fc:	jmp    94e5d8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x828>
  94e5fe:	mov    rdi,rbx
  94e601:	call   e1b80 <_Unwind_Resume@plt>

Disassembly of section .fini:
