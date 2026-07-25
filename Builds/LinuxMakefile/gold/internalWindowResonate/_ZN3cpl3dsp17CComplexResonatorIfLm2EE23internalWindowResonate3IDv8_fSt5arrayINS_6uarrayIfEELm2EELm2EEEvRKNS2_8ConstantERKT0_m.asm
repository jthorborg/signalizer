; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000953f30 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long):
  953f30:	endbr64
  953f34:	push   rbp
  953f35:	mov    rbp,rsp
  953f38:	push   r15
  953f3a:	push   r14
  953f3c:	mov    r14,rdi
  953f3f:	push   r13
  953f41:	push   r12
  953f43:	push   rbx
  953f44:	mov    rbx,rcx
  953f47:	and    rsp,0xffffffffffffffe0
  953f4b:	sub    rsp,0x180
  953f52:	mov    QWORD PTR [rsp+0x90],rsi
  953f5a:	mov    QWORD PTR [rsp+0x40],rdx
  953f5f:	movzx  r13d,BYTE PTR [rip+0x2df197]        # c330fe <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  953f67:	test   r13b,r13b
  953f6a:	je     954601 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x6d1>
  953f70:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  953f7c:	cmp    QWORD PTR [rax-0x1d8],0x0
  953f84:	lea    r15,[rax-0x360]
  953f8b:	je     953fb6 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x86>
  953f8d:	movzx  r12d,BYTE PTR [r15+0x180]
  953f95:	cmp    r12b,0xf
  953f99:	jbe    95466f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x73f>
  953f9f:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  953fab:	add    r12d,0x1
  953faf:	mov    BYTE PTR [rax-0x1e0],r12b
  953fb6:	mov    rax,QWORD PTR [rsp+0x90]
  953fbe:	mov    r9,QWORD PTR [r14+0x8]
  953fc2:	mov    rdi,QWORD PTR [r14]
  953fc5:	mov    rdx,QWORD PTR [rax+0x48]
  953fc9:	mov    rax,QWORD PTR [rax+0x38]
  953fcd:	mov    rcx,r9
  953fd0:	sub    rcx,rdi
  953fd3:	sar    rcx,0x2
  953fd7:	mov    r8,rax
  953fda:	imul   r8,rdx
  953fde:	lea    rsi,[r8*4+0x0]
  953fe6:	cmp    rcx,rsi
  953fe9:	jb     95464f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x71f>
  953fef:	cmp    rsi,rcx
  953ff2:	jb     9545e8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x6b8>
  953ff8:	mov    rdi,QWORD PTR [rsp+0x90]
  954000:	lea    rcx,[rdx+rdx*1]
  954004:	imul   rax,rcx
  954008:	cmp    QWORD PTR [rdi+0x40],0x0
  95400d:	je     954483 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x553>
  954013:	lea    r11,[rdx*4+0x0]
  95401b:	xor    r10d,r10d
  95401e:	mov    QWORD PTR [rsp+0x98],0x0
  95402a:	mov    r12,rbx
  95402d:	lea    rsi,[r11+rax*1]
  954031:	shl    rax,0x2
  954035:	mov    QWORD PTR [rsp+0x38],r11
  95403a:	mov    QWORD PTR [rsp+0x20],rax
  95403f:	lea    rdi,[rdx*8+0x0]
  954047:	mov    rax,rsi
  95404a:	shl    rdx,0x4
  95404e:	sub    rax,rcx
  954051:	mov    QWORD PTR [rsp+0x30],rdi
  954056:	shl    rax,0x2
  95405a:	mov    QWORD PTR [rsp+0x28],rdx
  95405f:	mov    QWORD PTR [rsp+0x18],rax
  954064:	lea    rax,[rsi*4+0x0]
  95406c:	mov    QWORD PTR [rsp+0x88],rax
  954074:	nop    DWORD PTR [rax+0x0]
  954078:	mov    rax,QWORD PTR [rsp+0x90]
  954080:	mov    rdi,QWORD PTR [rsp+0x30]
  954085:	mov    rbx,QWORD PTR [rsp+0x38]
  95408a:	mov    rcx,QWORD PTR [rsp+0x28]
  95408f:	mov    rsi,QWORD PTR [r14]
  954092:	mov    rax,QWORD PTR [rax]
  954095:	lea    r8,[r10+rdi*1]
  954099:	mov    rdx,QWORD PTR [rsp+0x40]
  95409e:	lea    rdi,[rbx+r8*1]
  9540a2:	lea    r9,[r11+rcx*1]
  9540a6:	add    rbx,rdi
  9540a9:	vmovups xmm1,XMMWORD PTR [rax+r8*1]
  9540af:	vinsertf128 ymm9,ymm1,XMMWORD PTR [rax+r8*1+0x10],0x1
  9540b7:	vmovups xmm2,XMMWORD PTR [rax+rdi*1]
  9540bc:	vinsertf128 ymm10,ymm2,XMMWORD PTR [rax+rdi*1+0x10],0x1
  9540c4:	vmovups xmm5,XMMWORD PTR [rax+rbx*1]
  9540c9:	vinsertf128 ymm11,ymm5,XMMWORD PTR [rax+rbx*1+0x10],0x1
  9540d1:	vmovups xmm1,XMMWORD PTR [rax+r9*1]
  9540d7:	vinsertf128 ymm8,ymm1,XMMWORD PTR [rax+r9*1+0x10],0x1
  9540df:	vmovups xmm6,XMMWORD PTR [rax+r10*1]
  9540e5:	vinsertf128 ymm6,ymm6,XMMWORD PTR [rax+r10*1+0x10],0x1
  9540ed:	vmovups xmm7,XMMWORD PTR [rax+r11*1]
  9540f3:	vinsertf128 ymm7,ymm7,XMMWORD PTR [rax+r11*1+0x10],0x1
  9540fb:	lea    rax,[rsi+r10*1]
  9540ff:	vmovups xmm2,XMMWORD PTR [rax]
  954103:	vinsertf128 ymm0,ymm2,XMMWORD PTR [rax+0x10],0x1
  95410a:	mov    rcx,QWORD PTR [rdx]
  95410d:	vmovups xmm2,XMMWORD PTR [rsi+r8*1]
  954113:	vinsertf128 ymm1,ymm2,XMMWORD PTR [rsi+r8*1+0x10],0x1
  95411b:	vmovaps YMMWORD PTR [rsp+0x120],ymm1
  954124:	vmovups xmm2,XMMWORD PTR [rsi+rdi*1]
  954129:	vinsertf128 ymm1,ymm2,XMMWORD PTR [rsi+rdi*1+0x10],0x1
  954131:	vmovaps YMMWORD PTR [rsp+0x100],ymm1
  95413a:	vmovups xmm2,XMMWORD PTR [rsi+rbx*1]
  95413f:	vinsertf128 ymm1,ymm2,XMMWORD PTR [rsi+rbx*1+0x10],0x1
  954147:	vmovups xmm5,XMMWORD PTR [rsi+r11*1]
  95414d:	vinsertf128 ymm5,ymm5,XMMWORD PTR [rsi+r11*1+0x10],0x1
  954155:	vmovaps YMMWORD PTR [rsp+0x140],ymm1
  95415e:	mov    rdx,QWORD PTR [rdx+0x10]
  954162:	vmovups xmm2,XMMWORD PTR [rsi+r9*1]
  954168:	vinsertf128 ymm1,ymm2,XMMWORD PTR [rsi+r9*1+0x10],0x1
  954170:	vmovaps YMMWORD PTR [rsp+0xc0],ymm1
  954179:	mov    QWORD PTR [rsp+0x48],rdx
  95417e:	mov    rdx,QWORD PTR [rsp+0x20]
  954183:	lea    r15,[r10+rdx*1]
  954187:	vmovups xmm1,XMMWORD PTR [rsi+r15*1]
  95418d:	mov    QWORD PTR [rsp+0x80],r15
  954195:	vinsertf128 ymm4,ymm1,XMMWORD PTR [rsi+r15*1+0x10],0x1
  95419d:	lea    r15,[r11+rdx*1]
  9541a1:	mov    rdx,QWORD PTR [rsp+0x18]
  9541a6:	mov    QWORD PTR [rsp+0x78],r15
  9541ab:	vmovups xmm1,XMMWORD PTR [rsi+r15*1]
  9541b1:	vinsertf128 ymm2,ymm1,XMMWORD PTR [rsi+r15*1+0x10],0x1
  9541b9:	vmovaps YMMWORD PTR [rsp+0xa0],ymm2
  9541c2:	lea    r15,[r10+rdx*1]
  9541c6:	lea    r13,[r11+rdx*1]
  9541ca:	mov    rdx,QWORD PTR [rsp+0x88]
  9541d2:	vmovups xmm2,XMMWORD PTR [rsi+r15*1]
  9541d8:	vinsertf128 ymm3,ymm2,XMMWORD PTR [rsi+r15*1+0x10],0x1
  9541e0:	vmovups xmm2,XMMWORD PTR [rsi+r13*1]
  9541e6:	vinsertf128 ymm1,ymm2,XMMWORD PTR [rsi+r13*1+0x10],0x1
  9541ee:	vmovaps YMMWORD PTR [rsp+0xe0],ymm1
  9541f7:	add    rdx,r10
  9541fa:	vmovups xmm1,XMMWORD PTR [rsi+rdx*1]
  9541ff:	mov    QWORD PTR [rsp+0x70],rdx
  954204:	vinsertf128 ymm2,ymm1,XMMWORD PTR [rsi+rdx*1+0x10],0x1
  95420c:	mov    rdx,QWORD PTR [rsp+0x88]
  954214:	vmovaps XMMWORD PTR [rsp+0x60],xmm1
  95421a:	add    rdx,r11
  95421d:	vmovups xmm1,XMMWORD PTR [rsi+rdx*1]
  954222:	mov    QWORD PTR [rsp+0x60],rdx
  954227:	vmovaps XMMWORD PTR [rsp+0x50],xmm1
  95422d:	vinsertf128 ymm1,ymm1,XMMWORD PTR [rsi+rdx*1+0x10],0x1
  954235:	test   r12,r12
  954238:	je     9543a7 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x477>
  95423e:	mov    rdx,QWORD PTR [rsp+0x48]
  954243:	mov    rax,rcx
  954246:	xor    ecx,ecx
  954248:	nop    DWORD PTR [rax+rax*1+0x0]
  954250:	vmulps ymm14,ymm0,ymm7
  954254:	vbroadcastss ymm12,DWORD PTR [rax]
  954259:	vmulps ymm0,ymm0,ymm6
  95425d:	vmulps ymm13,ymm7,ymm5
  954261:	vmulps ymm5,ymm6,ymm5
  954265:	vaddps ymm5,ymm14,ymm5
  954269:	add    rcx,0x1
  95426d:	add    rdx,0x4
  954271:	vaddps ymm0,ymm0,ymm12
  954276:	add    rax,0x4
  95427a:	vmovaps ymm14,YMMWORD PTR [rsp+0x100]
  954283:	vmulps ymm15,ymm10,YMMWORD PTR [rsp+0x120]
  95428c:	vsubps ymm0,ymm0,ymm13
  954291:	vmulps ymm13,ymm10,ymm14
  954296:	vmulps ymm14,ymm9,ymm14
  95429b:	vaddps ymm14,ymm15,ymm14
  9542a0:	vmulps ymm15,ymm8,YMMWORD PTR [rsp+0x140]
  9542a9:	vmovaps YMMWORD PTR [rsp+0x100],ymm14
  9542b2:	vmulps ymm14,ymm9,YMMWORD PTR [rsp+0x120]
  9542bb:	vaddps ymm14,ymm14,ymm12
  9542c0:	vsubps ymm13,ymm14,ymm13
  9542c5:	vmovaps ymm14,YMMWORD PTR [rsp+0xc0]
  9542ce:	vmovaps YMMWORD PTR [rsp+0x120],ymm13
  9542d7:	vmulps ymm13,ymm8,ymm14
  9542dc:	vmulps ymm14,ymm11,ymm14
  9542e1:	vaddps ymm14,ymm15,ymm14
  9542e6:	vmovaps ymm15,YMMWORD PTR [rsp+0xa0]
  9542ef:	vmovaps YMMWORD PTR [rsp+0xc0],ymm14
  9542f8:	vmulps ymm14,ymm11,YMMWORD PTR [rsp+0x140]
  954301:	vaddps ymm14,ymm14,ymm12
  954306:	vmulps ymm12,ymm6,ymm4
  95430a:	vmulps ymm4,ymm7,ymm4
  95430e:	vsubps ymm13,ymm14,ymm13
  954313:	vmulps ymm14,ymm7,ymm15
  954318:	vsubps ymm12,ymm12,ymm14
  95431d:	vmulps ymm14,ymm6,ymm15
  954322:	vmovaps ymm15,YMMWORD PTR [rsp+0xe0]
  95432b:	vaddps ymm4,ymm4,ymm14
  954330:	vmovaps YMMWORD PTR [rsp+0x140],ymm13
  954339:	vbroadcastss ymm13,DWORD PTR [rdx-0x4]
  95433f:	vmulps ymm14,ymm10,ymm15
  954344:	vmovaps YMMWORD PTR [rsp+0xa0],ymm4
  95434d:	vaddps ymm4,ymm13,ymm12
  954352:	vmulps ymm12,ymm3,ymm9
  954357:	vmulps ymm3,ymm3,ymm10
  95435c:	vsubps ymm12,ymm12,ymm14
  954361:	vmulps ymm14,ymm9,ymm15
  954366:	vaddps ymm3,ymm3,ymm14
  95436b:	vmulps ymm14,ymm2,ymm11
  954370:	vmulps ymm2,ymm8,ymm2
  954374:	vmovaps YMMWORD PTR [rsp+0xe0],ymm3
  95437d:	vaddps ymm3,ymm13,ymm12
  954382:	vmulps ymm12,ymm1,ymm8
  954387:	vmulps ymm1,ymm1,ymm11
  95438c:	vsubps ymm12,ymm14,ymm12
  954391:	vaddps ymm1,ymm2,ymm1
  954395:	vaddps ymm2,ymm13,ymm12
  95439a:	cmp    r12,rcx
  95439d:	jne    954250 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x320>
  9543a3:	lea    rax,[rsi+r10*1]
  9543a7:	vmovaps YMMWORD PTR [rax],ymm0
  9543ab:	mov    rax,QWORD PTR [r14]
  9543ae:	add    r10,0x20
  9543b2:	vmovaps ymm6,YMMWORD PTR [rsp+0x120]
  9543bb:	vmovaps ymm7,YMMWORD PTR [rsp+0x100]
  9543c4:	vmovaps YMMWORD PTR [rax+r11*1],ymm5
  9543ca:	mov    rax,QWORD PTR [r14]
  9543cd:	add    r11,0x20
  9543d1:	vmovaps YMMWORD PTR [rax+r8*1],ymm6
  9543d7:	mov    rax,QWORD PTR [r14]
  9543da:	vmovaps ymm6,YMMWORD PTR [rsp+0x140]
  9543e3:	vmovaps YMMWORD PTR [rax+rdi*1],ymm7
  9543e8:	mov    rax,QWORD PTR [r14]
  9543eb:	vmovaps ymm7,YMMWORD PTR [rsp+0xc0]
  9543f4:	mov    rdi,QWORD PTR [rsp+0x80]
  9543fc:	vmovaps YMMWORD PTR [rax+rbx*1],ymm6
  954401:	mov    rax,QWORD PTR [r14]
  954404:	vmovaps ymm6,YMMWORD PTR [rsp+0xa0]
  95440d:	mov    rbx,QWORD PTR [rsp+0x70]
  954412:	vmovaps YMMWORD PTR [rax+r9*1],ymm7
  954418:	mov    rax,QWORD PTR [r14]
  95441b:	vmovaps ymm7,YMMWORD PTR [rsp+0xe0]
  954424:	vmovaps YMMWORD PTR [rax+rdi*1],ymm4
  954429:	mov    rdi,QWORD PTR [rsp+0x78]
  95442e:	mov    rax,QWORD PTR [r14]
  954431:	vmovaps YMMWORD PTR [rax+rdi*1],ymm6
  954436:	mov    rax,QWORD PTR [r14]
  954439:	mov    rdi,QWORD PTR [rsp+0x60]
  95443e:	vmovaps YMMWORD PTR [rax+r15*1],ymm3
  954444:	mov    rax,QWORD PTR [r14]
  954447:	vmovaps YMMWORD PTR [rax+r13*1],ymm7
  95444d:	mov    rax,QWORD PTR [r14]
  954450:	vmovaps YMMWORD PTR [rax+rbx*1],ymm2
  954455:	mov    rax,QWORD PTR [r14]
  954458:	vmovaps YMMWORD PTR [rax+rdi*1],ymm1
  95445d:	mov    rbx,QWORD PTR [rsp+0x90]
  954465:	add    QWORD PTR [rsp+0x98],0x8
  95446e:	mov    rax,QWORD PTR [rsp+0x98]
  954476:	cmp    rax,QWORD PTR [rbx+0x40]
  95447a:	jb     954078 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x148>
  954480:	vzeroupper
  954483:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  95448f:	lea    r12,[rax-0x360]
  954496:	mov    r13,rax
  954499:	mov    r14,QWORD PTR [r12+0x188]
  9544a1:	test   r14,r14
  9544a4:	je     9544c8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x598>
  9544a6:	movzx  eax,BYTE PTR [r12+0x180]
  9544af:	lea    ebx,[rax-0x1]
  9544b2:	mov    BYTE PTR [r12+0x180],bl
  9544ba:	cmp    bl,0xf
  9544bd:	jbe    9544d7 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x5a7>
  9544bf:	add    WORD PTR [r14+0xfe2],0x1
  9544c8:	lea    rsp,[rbp-0x28]
  9544cc:	pop    rbx
  9544cd:	pop    r12
  9544cf:	pop    r13
  9544d1:	pop    r14
  9544d3:	pop    r15
  9544d5:	pop    rbp
  9544d6:	ret
  9544d7:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  9544dc:	movzx  r15d,bl
  9544e0:	movsxd rcx,r15d
  9544e3:	mov    rsi,rax
  9544e6:	vmovq  xmm6,rax
  9544eb:	lea    rdx,[rcx+rcx*2]
  9544ef:	shl    rdx,0x3
  9544f3:	sub    rsi,QWORD PTR [r12+rdx*1]
  9544f7:	vpinsrq xmm0,xmm6,rsi,0x1
  9544fd:	vpsubq xmm0,xmm0,XMMWORD PTR [r13+rdx*1-0x360]
  954507:	movzx  r13d,BYTE PTR [r12+0x181]
  954510:	cmp    r13b,bl
  954513:	jae    954529 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x5f9>
  954515:	lea    eax,[r15-0x1]
  954519:	vmovq  rdx,xmm0
  95451e:	cdqe
  954520:	lea    rax,[rax+rax*2]
  954524:	add    QWORD PTR [r12+rax*8+0x8],rdx
  954529:	movzx  r12d,WORD PTR [r14+0xfe0]
  954531:	cmp    r12w,0x7f
  954536:	je     9546cc <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x79c>
  95453c:	mov    QWORD PTR [rsp+0x120],rcx
  954544:	sub    ebx,r13d
  954547:	vmovdqa XMMWORD PTR [rsp+0x140],xmm0
  954550:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  95455c:	mov    rcx,QWORD PTR [rsp+0x120]
  954564:	mov    DWORD PTR [rsp+0x178],0x0
  95456f:	vmovdqa xmm0,XMMWORD PTR [rsp+0x140]
  954578:	vmovdqu XMMWORD PTR [rsp+0x168],xmm0
  954581:	lea    rdx,[rcx+rcx*2]
  954585:	lea    rdx,[rax+rdx*8-0x360]
  95458d:	mov    rax,QWORD PTR [rdx]
  954590:	mov    QWORD PTR [rsp+0x160],rax
  954598:	movzx  eax,BYTE PTR [rdx+0x10]
  95459c:	lea    edx,[r12+0x1]
  9545a1:	vmovdqa xmm7,XMMWORD PTR [rsp+0x160]
  9545aa:	mov    WORD PTR [r14+0xfe0],dx
  9545b2:	mov    ah,bl
  9545b4:	mov    WORD PTR [rsp+0x17c],ax
  9545bc:	movzx  eax,r12w
  9545c0:	vmovdqu xmm6,XMMWORD PTR [rsp+0x16e]
  9545c9:	shl    rax,0x5
  9545cd:	add    rax,r14
  9545d0:	vmovdqu XMMWORD PTR [rax],xmm7
  9545d4:	vmovdqu XMMWORD PTR [rax+0xe],xmm6
  9545d9:	lea    rsp,[rbp-0x28]
  9545dd:	pop    rbx
  9545de:	pop    r12
  9545e0:	pop    r13
  9545e2:	pop    r14
  9545e4:	pop    r15
  9545e6:	pop    rbp
  9545e7:	ret
  9545e8:	shl    r8,0x4
  9545ec:	add    rdi,r8
  9545ef:	cmp    r9,rdi
  9545f2:	je     953ff8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc8>
  9545f8:	mov    QWORD PTR [r14+0x8],rdi
  9545fc:	jmp    953ff8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc8>
  954601:	mov    eax,0x1
  954606:	lock xadd DWORD PTR [rip+0x2d3552],eax        # c27b60 <cpl::Profiling::registerRegion(char const*)::counter>
  95460e:	mov    edx,0x1
  954613:	mov    r13d,0x1
  954619:	add    eax,0x2
  95461c:	cmp    eax,0xfe
  954621:	ja     95463c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x70c>
  954623:	lea    rdx,[rip+0x2d19f6]        # c26020 <cpl::Profiling::regions>
  95462a:	mov    ecx,eax
  95462c:	mov    r13d,eax
  95462f:	lea    rdi,[rip+0xacc7c]        # a012b2 <_IO_stdin_used+0x182b2>
  954636:	mov    QWORD PTR [rdx+rcx*8],rdi
  95463a:	mov    edx,eax
  95463c:	xor    eax,eax
  95463e:	lock cmpxchg BYTE PTR [rip+0x2deab8],dl        # c330fe <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  954646:	cmovne r13d,eax
  95464a:	jmp    953f70 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x40>
  95464f:	sub    rsi,rcx
  954652:	mov    rdi,r14
  954655:	call   939bb0 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  95465a:	mov    rdi,QWORD PTR [rsp+0x90]
  954662:	mov    rax,QWORD PTR [rdi+0x38]
  954666:	mov    rdx,QWORD PTR [rdi+0x48]
  95466a:	jmp    953ff8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc8>
  95466f:	mov    QWORD PTR [rsp+0x140],rax
  954677:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  95467c:	mov    rdx,QWORD PTR [rsp+0x140]
  954684:	mov    QWORD PTR [rsp+0x168],0x0
  954690:	mov    QWORD PTR [rsp+0x160],rax
  954698:	lea    rax,[r12+r12*2]
  95469c:	vmovdqa xmm7,XMMWORD PTR [rsp+0x160]
  9546a5:	shl    rax,0x3
  9546a9:	mov    BYTE PTR [rsp+0x170],r13b
  9546b1:	vmovdqu XMMWORD PTR [rdx+rax*1-0x360],xmm7
  9546ba:	mov    BYTE PTR [rax+r15*1+0x10],r13b
  9546bf:	movzx  r12d,BYTE PTR [r15+0x180]
  9546c7:	jmp    953f9f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x6f>
  9546cc:	movzx  eax,WORD PTR [r14+0xfe2]
  9546d4:	cmp    ax,0xffff
  9546d8:	je     9544c8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x598>
  9546de:	add    eax,0x1
  9546e1:	mov    WORD PTR [r14+0xfe2],ax
  9546e9:	jmp    9544c8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x598>
  9546ee:	endbr64
  9546f2:	mov    rbx,rax
  9546f5:	vzeroupper
  9546f8:	call   93c790 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  9546fd:	mov    rdi,rbx
  954700:	call   e1b10 <_Unwind_Resume@plt>

Disassembly of section .fini:
