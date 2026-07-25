; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000944380 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long):
  944380:	endbr64
  944384:	push   r15
  944386:	push   r14
  944388:	push   r13
  94438a:	mov    r13,rdi
  94438d:	push   r12
  94438f:	push   rbp
  944390:	push   rbx
  944391:	mov    rbx,rcx
  944394:	sub    rsp,0x118
  94439b:	mov    QWORD PTR [rsp+0x88],rsi
  9443a3:	mov    QWORD PTR [rsp+0xe8],rdx
  9443ab:	movzx  r14d,BYTE PTR [rip+0x2dad3a]        # c1f0ed <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  9443b3:	test   r14b,r14b
  9443b6:	je     94489f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x51f>
  9443bc:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  9443c8:	cmp    QWORD PTR [rax-0x1d8],0x0
  9443d0:	mov    r12,rax
  9443d3:	lea    r15,[rax-0x360]
  9443da:	je     944404 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x84>
  9443dc:	movzx  ebp,BYTE PTR [r15+0x180]
  9443e4:	cmp    bpl,0xf
  9443e8:	jbe    944909 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x589>
  9443ee:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  9443fa:	add    ebp,0x1
  9443fd:	mov    BYTE PTR [rax-0x1e0],bpl
  944404:	mov    rax,QWORD PTR [rsp+0x88]
  94440c:	mov    r8,QWORD PTR [r13+0x8]
  944410:	mov    rcx,QWORD PTR [r13+0x0]
  944414:	mov    rdx,QWORD PTR [rax+0x48]
  944418:	mov    rdi,QWORD PTR [rax+0x38]
  94441c:	mov    rax,r8
  94441f:	sub    rax,rcx
  944422:	imul   rdi,rdx
  944426:	sar    rax,0x2
  94442a:	lea    rsi,[rdi*4+0x0]
  944432:	cmp    rax,rsi
  944435:	jb     9448ed <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x56d>
  94443b:	cmp    rsi,rax
  94443e:	jb     944886 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x506>
  944444:	mov    rax,QWORD PTR [rsp+0x88]
  94444c:	cmp    QWORD PTR [rax+0x40],0x0
  944451:	je     944740 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3c0>
  944457:	lea    r15,[rdx*4+0x0]
  94445f:	lea    r14,[rdx+rdx*2]
  944463:	mov    QWORD PTR [rsp+0x60],0x0
  94446c:	mov    QWORD PTR [rsp+0x68],r15
  944471:	shl    r14,0x2
  944475:	nop    DWORD PTR [rax]
  944478:	mov    rax,QWORD PTR [rsp+0x88]
  944480:	mov    rbp,QWORD PTR [rax]
  944483:	mov    rax,QWORD PTR [rsp+0x60]
  944488:	shl    rax,0x2
  94448c:	lea    rdi,[rbp+rax*1+0x0]
  944491:	mov    QWORD PTR [rsp+0xe0],rax
  944499:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  94449e:	lea    rdi,[rbp+r15*1+0x0]
  9444a3:	movaps XMMWORD PTR [rsp],xmm0
  9444a7:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  9444ac:	mov    r12,QWORD PTR [rsp+0x68]
  9444b1:	movaps XMMWORD PTR [rsp+0x10],xmm0
  9444b6:	lea    rax,[r12+r15*1]
  9444ba:	lea    r12,[r12+r14*1]
  9444be:	lea    rdi,[rbp+rax*1+0x0]
  9444c3:	mov    QWORD PTR [rsp+0x78],rax
  9444c8:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  9444cd:	lea    rdi,[rbp+r14*1+0x0]
  9444d2:	movaps XMMWORD PTR [rsp+0x20],xmm0
  9444d7:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  9444dc:	lea    rdi,[rbp+r12*1+0x0]
  9444e1:	movaps XMMWORD PTR [rsp+0x30],xmm0
  9444e6:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  9444eb:	mov    rcx,QWORD PTR [rsp+0x68]
  9444f0:	mov    QWORD PTR [rsp+0x80],r12
  9444f8:	movaps XMMWORD PTR [rsp+0x40],xmm0
  9444fd:	lea    rdx,[rcx+r12*1]
  944501:	lea    rdi,[rbp+rdx*1+0x0]
  944506:	mov    QWORD PTR [rsp+0x70],rdx
  94450b:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944510:	mov    rax,QWORD PTR [rsp+0xe0]
  944518:	mov    r12,QWORD PTR [r13+0x0]
  94451c:	mov    rcx,QWORD PTR [rsp+0xe8]
  944524:	movaps XMMWORD PTR [rsp+0x50],xmm0
  944529:	mov    rbp,QWORD PTR [rcx]
  94452c:	lea    rcx,[r12+rax*1]
  944530:	mov    rdi,rcx
  944533:	mov    QWORD PTR [rsp+0xe0],rcx
  94453b:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944540:	lea    rdi,[r12+r15*1]
  944544:	movaps XMMWORD PTR [rsp+0xd0],xmm0
  94454c:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944551:	mov    rsi,QWORD PTR [rsp+0x78]
  944556:	movaps XMMWORD PTR [rsp+0xc0],xmm0
  94455e:	lea    rdi,[r12+rsi*1]
  944562:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944567:	lea    rdi,[r12+r14*1]
  94456b:	movaps XMMWORD PTR [rsp+0xb0],xmm0
  944573:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944578:	mov    rax,QWORD PTR [rsp+0x80]
  944580:	movaps XMMWORD PTR [rsp+0xa0],xmm0
  944588:	lea    rdi,[r12+rax*1]
  94458c:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  944591:	mov    rdx,QWORD PTR [rsp+0x70]
  944596:	movaps XMMWORD PTR [rsp+0x90],xmm0
  94459e:	lea    rdi,[r12+rdx*1]
  9445a2:	call   85faf0 <float __vector(4) cpl::simd::load<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*)>
  9445a7:	test   rbx,rbx
  9445aa:	movaps xmm2,XMMWORD PTR [rsp+0x90]
  9445b2:	movaps xmm4,XMMWORD PTR [rsp+0xa0]
  9445ba:	movaps xmm5,XMMWORD PTR [rsp+0xb0]
  9445c2:	movaps xmm6,XMMWORD PTR [rsp+0xc0]
  9445ca:	movaps xmm1,xmm0
  9445cd:	movaps xmm7,XMMWORD PTR [rsp+0xd0]
  9445d5:	je     944687 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x307>
  9445db:	xor    eax,eax
  9445dd:	nop    DWORD PTR [rax]
  9445e0:	movss  xmm0,DWORD PTR [rbp+rax*4+0x0]
  9445e6:	add    rax,0x1
  9445ea:	call   9205e0 <float __vector(4) cpl::simd::broadcast<float __vector(4)>(cpl::simd::scalar_of<float __vector(4), 4ul>::type const*) [clone .isra.0]>
  9445ef:	movaps xmm9,XMMWORD PTR [rsp]
  9445f4:	movaps xmm10,XMMWORD PTR [rsp+0x10]
  9445fa:	movaps xmm11,XMMWORD PTR [rsp+0x20]
  944600:	movaps xmm12,XMMWORD PTR [rsp+0x30]
  944606:	movaps xmm3,xmm9
  94460a:	movaps xmm8,xmm10
  94460e:	movaps xmm13,XMMWORD PTR [rsp+0x40]
  944614:	movaps xmm14,XMMWORD PTR [rsp+0x50]
  94461a:	mulps  xmm8,xmm6
  94461e:	mulps  xmm3,xmm7
  944621:	mulps  xmm6,xmm9
  944625:	mulps  xmm7,xmm10
  944629:	subps  xmm3,xmm8
  94462d:	movaps xmm8,xmm12
  944631:	mulps  xmm8,xmm4
  944635:	addps  xmm6,xmm7
  944638:	mulps  xmm4,xmm11
  94463c:	movaps xmm7,xmm3
  94463f:	movaps xmm3,xmm11
  944643:	mulps  xmm3,xmm5
  944646:	addps  xmm7,xmm0
  944649:	mulps  xmm5,xmm12
  94464d:	subps  xmm3,xmm8
  944651:	movaps xmm8,xmm14
  944655:	addps  xmm4,xmm5
  944658:	mulps  xmm8,xmm1
  94465c:	mulps  xmm1,xmm13
  944660:	movaps xmm5,xmm3
  944663:	movaps xmm3,xmm13
  944667:	mulps  xmm3,xmm2
  94466a:	addps  xmm5,xmm0
  94466d:	mulps  xmm2,xmm14
  944671:	subps  xmm3,xmm8
  944675:	addps  xmm1,xmm2
  944678:	movaps xmm2,xmm0
  94467b:	addps  xmm2,xmm3
  94467e:	cmp    rbx,rax
  944681:	jne    9445e0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x260>
  944687:	mov    rdi,QWORD PTR [rsp+0xe0]
  94468f:	movaps xmm0,xmm7
  944692:	movaps XMMWORD PTR [rsp],xmm1
  944696:	movaps XMMWORD PTR [rsp+0x10],xmm2
  94469b:	movaps XMMWORD PTR [rsp+0x30],xmm4
  9446a0:	movaps XMMWORD PTR [rsp+0x40],xmm5
  9446a5:	movaps XMMWORD PTR [rsp+0x20],xmm6
  9446aa:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  9446af:	movaps xmm6,XMMWORD PTR [rsp+0x20]
  9446b4:	mov    rdi,QWORD PTR [r13+0x0]
  9446b8:	movaps xmm0,xmm6
  9446bb:	add    rdi,r15
  9446be:	add    r15,0x10
  9446c2:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  9446c7:	movaps xmm5,XMMWORD PTR [rsp+0x40]
  9446cc:	mov    rdi,QWORD PTR [rsp+0x78]
  9446d1:	add    rdi,QWORD PTR [r13+0x0]
  9446d5:	movaps xmm0,xmm5
  9446d8:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  9446dd:	movaps xmm4,XMMWORD PTR [rsp+0x30]
  9446e2:	mov    rdi,QWORD PTR [r13+0x0]
  9446e6:	movaps xmm0,xmm4
  9446e9:	add    rdi,r14
  9446ec:	add    r14,0x10
  9446f0:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  9446f5:	movaps xmm2,XMMWORD PTR [rsp+0x10]
  9446fa:	mov    rdi,QWORD PTR [rsp+0x80]
  944702:	add    rdi,QWORD PTR [r13+0x0]
  944706:	movaps xmm0,xmm2
  944709:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  94470e:	movaps xmm1,XMMWORD PTR [rsp]
  944712:	mov    rdi,QWORD PTR [rsp+0x70]
  944717:	add    rdi,QWORD PTR [r13+0x0]
  94471b:	movaps xmm0,xmm1
  94471e:	call   85fb20 <cpl::simd::store(float*, float __vector(4))>
  944723:	mov    rsi,QWORD PTR [rsp+0x88]
  94472b:	add    QWORD PTR [rsp+0x60],0x4
  944731:	mov    rax,QWORD PTR [rsp+0x60]
  944736:	cmp    rax,QWORD PTR [rsi+0x40]
  94473a:	jb     944478 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xf8>
  944740:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94474c:	lea    rbp,[rax-0x360]
  944753:	mov    r12,rax
  944756:	mov    r13,QWORD PTR [rbp+0x188]
  94475d:	test   r13,r13
  944760:	je     944780 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x400>
  944762:	movzx  eax,BYTE PTR [rbp+0x180]
  944769:	lea    ebx,[rax-0x1]
  94476c:	mov    BYTE PTR [rbp+0x180],bl
  944772:	cmp    bl,0xf
  944775:	jbe    944792 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x412>
  944777:	add    WORD PTR [r13+0xfe2],0x1
  944780:	add    rsp,0x118
  944787:	pop    rbx
  944788:	pop    rbp
  944789:	pop    r12
  94478b:	pop    r13
  94478d:	pop    r14
  94478f:	pop    r15
  944791:	ret
  944792:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  944797:	movzx  r14d,bl
  94479b:	lea    rdx,[r14+r14*2]
  94479f:	mov    rcx,rax
  9447a2:	movq   xmm0,rax
  9447a7:	shl    rdx,0x3
  9447ab:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  9447b0:	movq   xmm1,rcx
  9447b5:	punpcklqdq xmm0,xmm1
  9447b9:	movdqu xmm1,XMMWORD PTR [r12+rdx*1-0x360]
  9447c3:	movzx  r12d,BYTE PTR [rbp+0x181]
  9447cb:	psubq  xmm0,xmm1
  9447cf:	cmp    r12b,bl
  9447d2:	jae    9447f2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x472>
  9447d4:	lea    eax,[r14-0x1]
  9447d8:	movdqa xmm1,xmm0
  9447dc:	cdqe
  9447de:	lea    rax,[rax+rax*2]
  9447e2:	movq   xmm2,QWORD PTR [rbp+rax*8+0x8]
  9447e8:	paddq  xmm1,xmm2
  9447ec:	movq   QWORD PTR [rbp+rax*8+0x8],xmm1
  9447f2:	movzx  ebp,WORD PTR [r13+0xfe0]
  9447fa:	cmp    bp,0x7f
  9447fe:	je     944957 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x5d7>
  944804:	movaps XMMWORD PTR [rsp],xmm0
  944808:	sub    ebx,r12d
  94480b:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  944817:	lea    rdx,[r14+r14*2]
  94481b:	movdqa xmm0,XMMWORD PTR [rsp]
  944820:	mov    DWORD PTR [rsp+0x108],0x0
  94482b:	movups XMMWORD PTR [rsp+0xf8],xmm0
  944833:	lea    rdx,[rax+rdx*8-0x360]
  94483b:	mov    rax,QWORD PTR [rdx]
  94483e:	mov    QWORD PTR [rsp+0xf0],rax
  944846:	movzx  eax,BYTE PTR [rdx+0x10]
  94484a:	lea    edx,[rbp+0x1]
  94484d:	shl    rbp,0x5
  944851:	movdqa xmm1,XMMWORD PTR [rsp+0xf0]
  94485a:	mov    WORD PTR [r13+0xfe0],dx
  944862:	mov    ah,bl
  944864:	mov    WORD PTR [rsp+0x10c],ax
  94486c:	lea    rax,[r13+rbp*1+0x0]
  944871:	movups XMMWORD PTR [rax],xmm1
  944874:	movdqu xmm1,XMMWORD PTR [rsp+0xfe]
  94487d:	movups XMMWORD PTR [rax+0xe],xmm1
  944881:	jmp    944780 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x400>
  944886:	shl    rdi,0x4
  94488a:	add    rcx,rdi
  94488d:	cmp    r8,rcx
  944890:	je     944444 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc4>
  944896:	mov    QWORD PTR [r13+0x8],rcx
  94489a:	jmp    944444 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc4>
  94489f:	mov    eax,0x1
  9448a4:	lock xadd DWORD PTR [rip+0x2cf2b4],eax        # c13b60 <cpl::Profiling::registerRegion(char const*)::counter>
  9448ac:	mov    edx,0x1
  9448b1:	mov    r14d,0x1
  9448b7:	add    eax,0x2
  9448ba:	cmp    eax,0xfe
  9448bf:	ja     9448da <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x55a>
  9448c1:	lea    rdx,[rip+0x2cd758]        # c12020 <cpl::Profiling::regions>
  9448c8:	mov    ecx,eax
  9448ca:	lea    rsi,[rip+0xa3875]        # 9e8146 <_IO_stdin_used+0x18146>
  9448d1:	mov    r14d,eax
  9448d4:	mov    QWORD PTR [rdx+rcx*8],rsi
  9448d8:	mov    edx,eax
  9448da:	xor    eax,eax
  9448dc:	lock cmpxchg BYTE PTR [rip+0x2da809],dl        # c1f0ed <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  9448e4:	cmovne r14d,eax
  9448e8:	jmp    9443bc <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3c>
  9448ed:	sub    rsi,rax
  9448f0:	mov    rdi,r13
  9448f3:	call   91d820 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  9448f8:	mov    rax,QWORD PTR [rsp+0x88]
  944900:	mov    rdx,QWORD PTR [rax+0x48]
  944904:	jmp    944444 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc4>
  944909:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  94490e:	mov    BYTE PTR [rsp+0x100],r14b
  944916:	mov    QWORD PTR [rsp+0xf0],rax
  94491e:	lea    rax,[rbp+rbp*2+0x0]
  944923:	mov    QWORD PTR [rsp+0xf8],0x0
  94492f:	shl    rax,0x3
  944933:	movdqa xmm1,XMMWORD PTR [rsp+0xf0]
  94493c:	movups XMMWORD PTR [r12+rax*1-0x360],xmm1
  944945:	mov    BYTE PTR [rax+r15*1+0x10],r14b
  94494a:	movzx  ebp,BYTE PTR [r15+0x180]
  944952:	jmp    9443ee <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x6e>
  944957:	movzx  eax,WORD PTR [r13+0xfe2]
  94495f:	cmp    ax,0xffff
  944963:	je     944780 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x400>
  944969:	add    eax,0x1
  94496c:	mov    WORD PTR [r13+0xfe2],ax
  944974:	jmp    944780 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x400>
  944979:	endbr64
  94497d:	mov    rbx,rax
  944980:	call   920350 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  944985:	mov    rdi,rbx
  944988:	call   e1b80 <_Unwind_Resume@plt>

Disassembly of section .fini:
