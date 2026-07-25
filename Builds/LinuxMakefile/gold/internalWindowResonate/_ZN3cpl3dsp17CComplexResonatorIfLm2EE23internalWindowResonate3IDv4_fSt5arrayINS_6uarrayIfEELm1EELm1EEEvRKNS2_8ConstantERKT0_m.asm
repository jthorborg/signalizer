; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

000000000093f1a0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long):
  93f1a0:	endbr64
  93f1a4:	push   r15
  93f1a6:	push   r14
  93f1a8:	push   r13
  93f1aa:	mov    r13,rdi
  93f1ad:	push   r12
  93f1af:	push   rbp
  93f1b0:	push   rbx
  93f1b1:	mov    rbx,rcx
  93f1b4:	sub    rsp,0x118
  93f1bb:	mov    QWORD PTR [rsp+0x88],rsi
  93f1c3:	mov    QWORD PTR [rsp+0xe8],rdx
  93f1cb:	movzx  r14d,BYTE PTR [rip+0x2dff24]        # c1f0f7 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  93f1d3:	test   r14b,r14b
  93f1d6:	je     93f6bf <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x51f>
  93f1dc:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  93f1e8:	cmp    QWORD PTR [rax-0x1d8],0x0
  93f1f0:	mov    r12,rax
  93f1f3:	lea    r15,[rax-0x360]
  93f1fa:	je     93f224 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x84>
  93f1fc:	movzx  ebp,BYTE PTR [r15+0x180]
  93f204:	cmp    bpl,0xf
  93f208:	jbe    93f729 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x589>
  93f20e:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  93f21a:	add    ebp,0x1
  93f21d:	mov    BYTE PTR [rax-0x1e0],bpl
  93f224:	mov    rax,QWORD PTR [rsp+0x88]
  93f22c:	mov    r8,QWORD PTR [r13+0x8]
  93f230:	mov    rcx,QWORD PTR [r13+0x0]
  93f234:	mov    rdx,QWORD PTR [rax+0x48]
  93f238:	mov    rdi,QWORD PTR [rax+0x38]
  93f23c:	mov    rax,r8
  93f23f:	sub    rax,rcx
  93f242:	imul   rdi,rdx
  93f246:	sar    rax,0x2
  93f24a:	lea    rsi,[rdi*4+0x0]
  93f252:	cmp    rax,rsi
  93f255:	jb     93f70d <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x56d>
  93f25b:	cmp    rsi,rax
  93f25e:	jb     93f6a6 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x506>
  93f264:	mov    rax,QWORD PTR [rsp+0x88]
  93f26c:	cmp    QWORD PTR [rax+0x40],0x0
  93f271:	je     93f560 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3c0>
  93f277:	lea    r15,[rdx*4+0x0]
  93f27f:	lea    r14,[rdx+rdx*2]
  93f283:	mov    QWORD PTR [rsp+0x60],0x0
  93f28c:	mov    QWORD PTR [rsp+0x68],r15
  93f291:	shl    r14,0x2
  93f295:	nop    DWORD PTR [rax]
  93f298:	mov    rax,QWORD PTR [rsp+0x88]
  93f2a0:	mov    rbp,QWORD PTR [rax]
  93f2a3:	mov    rax,QWORD PTR [rsp+0x60]
  93f2a8:	shl    rax,0x2
  93f2ac:	lea    rdi,[rbp+rax*1+0x0]
  93f2b1:	mov    QWORD PTR [rsp+0xe0],rax
  93f2b9:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f2be:	lea    rdi,[rbp+r15*1+0x0]
  93f2c3:	movaps XMMWORD PTR [rsp],xmm0
  93f2c7:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f2cc:	mov    r12,QWORD PTR [rsp+0x68]
  93f2d1:	movaps XMMWORD PTR [rsp+0x10],xmm0
  93f2d6:	lea    rax,[r12+r15*1]
  93f2da:	lea    r12,[r12+r14*1]
  93f2de:	lea    rdi,[rbp+rax*1+0x0]
  93f2e3:	mov    QWORD PTR [rsp+0x78],rax
  93f2e8:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f2ed:	lea    rdi,[rbp+r14*1+0x0]
  93f2f2:	movaps XMMWORD PTR [rsp+0x20],xmm0
  93f2f7:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f2fc:	lea    rdi,[rbp+r12*1+0x0]
  93f301:	movaps XMMWORD PTR [rsp+0x30],xmm0
  93f306:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f30b:	mov    rcx,QWORD PTR [rsp+0x68]
  93f310:	mov    QWORD PTR [rsp+0x80],r12
  93f318:	movaps XMMWORD PTR [rsp+0x40],xmm0
  93f31d:	lea    rdx,[rcx+r12*1]
  93f321:	lea    rdi,[rbp+rdx*1+0x0]
  93f326:	mov    QWORD PTR [rsp+0x70],rdx
  93f32b:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f330:	mov    rax,QWORD PTR [rsp+0xe0]
  93f338:	mov    r12,QWORD PTR [r13+0x0]
  93f33c:	mov    rcx,QWORD PTR [rsp+0xe8]
  93f344:	movaps XMMWORD PTR [rsp+0x50],xmm0
  93f349:	mov    rbp,QWORD PTR [rcx]
  93f34c:	lea    rcx,[r12+rax*1]
  93f350:	mov    rdi,rcx
  93f353:	mov    QWORD PTR [rsp+0xe0],rcx
  93f35b:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f360:	lea    rdi,[r12+r15*1]
  93f364:	movaps XMMWORD PTR [rsp+0xd0],xmm0
  93f36c:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f371:	mov    rsi,QWORD PTR [rsp+0x78]
  93f376:	movaps XMMWORD PTR [rsp+0xc0],xmm0
  93f37e:	lea    rdi,[r12+rsi*1]
  93f382:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f387:	lea    rdi,[r12+r14*1]
  93f38b:	movaps XMMWORD PTR [rsp+0xb0],xmm0
  93f393:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f398:	mov    rax,QWORD PTR [rsp+0x80]
  93f3a0:	movaps XMMWORD PTR [rsp+0xa0],xmm0
  93f3a8:	lea    rdi,[r12+rax*1]
  93f3ac:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f3b1:	mov    rdx,QWORD PTR [rsp+0x70]
  93f3b6:	movaps XMMWORD PTR [rsp+0x90],xmm0
  93f3be:	lea    rdi,[r12+rdx*1]
  93f3c2:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  93f3c7:	test   rbx,rbx
  93f3ca:	movaps xmm2,XMMWORD PTR [rsp+0x90]
  93f3d2:	movaps xmm4,XMMWORD PTR [rsp+0xa0]
  93f3da:	movaps xmm5,XMMWORD PTR [rsp+0xb0]
  93f3e2:	movaps xmm6,XMMWORD PTR [rsp+0xc0]
  93f3ea:	movaps xmm1,xmm0
  93f3ed:	movaps xmm7,XMMWORD PTR [rsp+0xd0]
  93f3f5:	je     93f4a7 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x307>
  93f3fb:	xor    eax,eax
  93f3fd:	nop    DWORD PTR [rax]
  93f400:	movss  xmm0,DWORD PTR [rbp+rax*4+0x0]
  93f406:	add    rax,0x1
  93f40a:	call   9205e0 <float __vector(4) cpl::simd::broadcast<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*) [clone .isra.0]>
  93f40f:	movaps xmm9,XMMWORD PTR [rsp]
  93f414:	movaps xmm10,XMMWORD PTR [rsp+0x10]
  93f41a:	movaps xmm11,XMMWORD PTR [rsp+0x20]
  93f420:	movaps xmm12,XMMWORD PTR [rsp+0x30]
  93f426:	movaps xmm3,xmm9
  93f42a:	movaps xmm8,xmm10
  93f42e:	movaps xmm13,XMMWORD PTR [rsp+0x40]
  93f434:	movaps xmm14,XMMWORD PTR [rsp+0x50]
  93f43a:	mulps  xmm8,xmm6
  93f43e:	mulps  xmm3,xmm7
  93f441:	mulps  xmm6,xmm9
  93f445:	mulps  xmm7,xmm10
  93f449:	subps  xmm3,xmm8
  93f44d:	movaps xmm8,xmm12
  93f451:	mulps  xmm8,xmm4
  93f455:	addps  xmm6,xmm7
  93f458:	mulps  xmm4,xmm11
  93f45c:	movaps xmm7,xmm3
  93f45f:	movaps xmm3,xmm11
  93f463:	mulps  xmm3,xmm5
  93f466:	addps  xmm7,xmm0
  93f469:	mulps  xmm5,xmm12
  93f46d:	subps  xmm3,xmm8
  93f471:	movaps xmm8,xmm14
  93f475:	addps  xmm4,xmm5
  93f478:	mulps  xmm8,xmm1
  93f47c:	mulps  xmm1,xmm13
  93f480:	movaps xmm5,xmm3
  93f483:	movaps xmm3,xmm13
  93f487:	mulps  xmm3,xmm2
  93f48a:	addps  xmm5,xmm0
  93f48d:	mulps  xmm2,xmm14
  93f491:	subps  xmm3,xmm8
  93f495:	addps  xmm1,xmm2
  93f498:	movaps xmm2,xmm0
  93f49b:	addps  xmm2,xmm3
  93f49e:	cmp    rbx,rax
  93f4a1:	jne    93f400 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x260>
  93f4a7:	mov    rdi,QWORD PTR [rsp+0xe0]
  93f4af:	movaps xmm0,xmm7
  93f4b2:	movaps XMMWORD PTR [rsp],xmm1
  93f4b6:	movaps XMMWORD PTR [rsp+0x10],xmm2
  93f4bb:	movaps XMMWORD PTR [rsp+0x30],xmm4
  93f4c0:	movaps XMMWORD PTR [rsp+0x40],xmm5
  93f4c5:	movaps XMMWORD PTR [rsp+0x20],xmm6
  93f4ca:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  93f4cf:	movaps xmm6,XMMWORD PTR [rsp+0x20]
  93f4d4:	mov    rdi,QWORD PTR [r13+0x0]
  93f4d8:	movaps xmm0,xmm6
  93f4db:	add    rdi,r15
  93f4de:	add    r15,0x10
  93f4e2:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  93f4e7:	movaps xmm5,XMMWORD PTR [rsp+0x40]
  93f4ec:	mov    rdi,QWORD PTR [rsp+0x78]
  93f4f1:	add    rdi,QWORD PTR [r13+0x0]
  93f4f5:	movaps xmm0,xmm5
  93f4f8:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  93f4fd:	movaps xmm4,XMMWORD PTR [rsp+0x30]
  93f502:	mov    rdi,QWORD PTR [r13+0x0]
  93f506:	movaps xmm0,xmm4
  93f509:	add    rdi,r14
  93f50c:	add    r14,0x10
  93f510:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  93f515:	movaps xmm2,XMMWORD PTR [rsp+0x10]
  93f51a:	mov    rdi,QWORD PTR [rsp+0x80]
  93f522:	add    rdi,QWORD PTR [r13+0x0]
  93f526:	movaps xmm0,xmm2
  93f529:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  93f52e:	movaps xmm1,XMMWORD PTR [rsp]
  93f532:	mov    rdi,QWORD PTR [rsp+0x70]
  93f537:	add    rdi,QWORD PTR [r13+0x0]
  93f53b:	movaps xmm0,xmm1
  93f53e:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  93f543:	mov    rsi,QWORD PTR [rsp+0x88]
  93f54b:	add    QWORD PTR [rsp+0x60],0x4
  93f551:	mov    rax,QWORD PTR [rsp+0x60]
  93f556:	cmp    rax,QWORD PTR [rsi+0x40]
  93f55a:	jb     93f298 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xf8>
  93f560:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  93f56c:	lea    rbp,[rax-0x360]
  93f573:	mov    r12,rax
  93f576:	mov    r13,QWORD PTR [rbp+0x188]
  93f57d:	test   r13,r13
  93f580:	je     93f5a0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x400>
  93f582:	movzx  eax,BYTE PTR [rbp+0x180]
  93f589:	lea    ebx,[rax-0x1]
  93f58c:	mov    BYTE PTR [rbp+0x180],bl
  93f592:	cmp    bl,0xf
  93f595:	jbe    93f5b2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x412>
  93f597:	add    WORD PTR [r13+0xfe2],0x1
  93f5a0:	add    rsp,0x118
  93f5a7:	pop    rbx
  93f5a8:	pop    rbp
  93f5a9:	pop    r12
  93f5ab:	pop    r13
  93f5ad:	pop    r14
  93f5af:	pop    r15
  93f5b1:	ret
  93f5b2:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  93f5b7:	movzx  r14d,bl
  93f5bb:	lea    rdx,[r14+r14*2]
  93f5bf:	mov    rcx,rax
  93f5c2:	movq   xmm0,rax
  93f5c7:	shl    rdx,0x3
  93f5cb:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  93f5d0:	movq   xmm1,rcx
  93f5d5:	punpcklqdq xmm0,xmm1
  93f5d9:	movdqu xmm1,XMMWORD PTR [r12+rdx*1-0x360]
  93f5e3:	movzx  r12d,BYTE PTR [rbp+0x181]
  93f5eb:	psubq  xmm0,xmm1
  93f5ef:	cmp    r12b,bl
  93f5f2:	jae    93f612 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x472>
  93f5f4:	lea    eax,[r14-0x1]
  93f5f8:	movdqa xmm1,xmm0
  93f5fc:	cdqe
  93f5fe:	lea    rax,[rax+rax*2]
  93f602:	movq   xmm2,QWORD PTR [rbp+rax*8+0x8]
  93f608:	paddq  xmm1,xmm2
  93f60c:	movq   QWORD PTR [rbp+rax*8+0x8],xmm1
  93f612:	movzx  ebp,WORD PTR [r13+0xfe0]
  93f61a:	cmp    bp,0x7f
  93f61e:	je     93f777 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x5d7>
  93f624:	movaps XMMWORD PTR [rsp],xmm0
  93f628:	sub    ebx,r12d
  93f62b:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  93f637:	lea    rdx,[r14+r14*2]
  93f63b:	movdqa xmm0,XMMWORD PTR [rsp]
  93f640:	mov    DWORD PTR [rsp+0x108],0x0
  93f64b:	movups XMMWORD PTR [rsp+0xf8],xmm0
  93f653:	lea    rdx,[rax+rdx*8-0x360]
  93f65b:	mov    rax,QWORD PTR [rdx]
  93f65e:	mov    QWORD PTR [rsp+0xf0],rax
  93f666:	movzx  eax,BYTE PTR [rdx+0x10]
  93f66a:	lea    edx,[rbp+0x1]
  93f66d:	shl    rbp,0x5
  93f671:	movdqa xmm1,XMMWORD PTR [rsp+0xf0]
  93f67a:	mov    WORD PTR [r13+0xfe0],dx
  93f682:	mov    ah,bl
  93f684:	mov    WORD PTR [rsp+0x10c],ax
  93f68c:	lea    rax,[r13+rbp*1+0x0]
  93f691:	movups XMMWORD PTR [rax],xmm1
  93f694:	movdqu xmm1,XMMWORD PTR [rsp+0xfe]
  93f69d:	movups XMMWORD PTR [rax+0xe],xmm1
  93f6a1:	jmp    93f5a0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x400>
  93f6a6:	shl    rdi,0x4
  93f6aa:	add    rcx,rdi
  93f6ad:	cmp    r8,rcx
  93f6b0:	je     93f264 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc4>
  93f6b6:	mov    QWORD PTR [r13+0x8],rcx
  93f6ba:	jmp    93f264 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc4>
  93f6bf:	mov    eax,0x1
  93f6c4:	lock xadd DWORD PTR [rip+0x2d4494],eax        # c13b60 <cpl::Profiling::registerRegion(char const*)::counter>
  93f6cc:	mov    edx,0x1
  93f6d1:	mov    r14d,0x1
  93f6d7:	add    eax,0x2
  93f6da:	cmp    eax,0xfe
  93f6df:	ja     93f6fa <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x55a>
  93f6e1:	lea    rdx,[rip+0x2d2938]        # c12020 <cpl::Profiling::regions>
  93f6e8:	mov    ecx,eax
  93f6ea:	lea    rsi,[rip+0xa8a55]        # 9e8146 <_IO_stdin_used+0x18146>
  93f6f1:	mov    r14d,eax
  93f6f4:	mov    QWORD PTR [rdx+rcx*8],rsi
  93f6f8:	mov    edx,eax
  93f6fa:	xor    eax,eax
  93f6fc:	lock cmpxchg BYTE PTR [rip+0x2df9f3],dl        # c1f0f7 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  93f704:	cmovne r14d,eax
  93f708:	jmp    93f1dc <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3c>
  93f70d:	sub    rsi,rax
  93f710:	mov    rdi,r13
  93f713:	call   91d820 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  93f718:	mov    rax,QWORD PTR [rsp+0x88]
  93f720:	mov    rdx,QWORD PTR [rax+0x48]
  93f724:	jmp    93f264 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc4>
  93f729:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  93f72e:	mov    BYTE PTR [rsp+0x100],r14b
  93f736:	mov    QWORD PTR [rsp+0xf0],rax
  93f73e:	lea    rax,[rbp+rbp*2+0x0]
  93f743:	mov    QWORD PTR [rsp+0xf8],0x0
  93f74f:	shl    rax,0x3
  93f753:	movdqa xmm1,XMMWORD PTR [rsp+0xf0]
  93f75c:	movups XMMWORD PTR [r12+rax*1-0x360],xmm1
  93f765:	mov    BYTE PTR [rax+r15*1+0x10],r14b
  93f76a:	movzx  ebp,BYTE PTR [r15+0x180]
  93f772:	jmp    93f20e <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x6e>
  93f777:	movzx  eax,WORD PTR [r13+0xfe2]
  93f77f:	cmp    ax,0xffff
  93f783:	je     93f5a0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x400>
  93f789:	add    eax,0x1
  93f78c:	mov    WORD PTR [r13+0xfe2],ax
  93f794:	jmp    93f5a0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x400>
  93f799:	endbr64
  93f79d:	mov    rbx,rax
  93f7a0:	call   920350 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  93f7a5:	mov    rdi,rbx
  93f7a8:	call   e1b80 <_Unwind_Resume@plt>

Disassembly of section .fini:
