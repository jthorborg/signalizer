; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000949560 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long):
  949560:	endbr64
  949564:	push   r15
  949566:	push   r14
  949568:	push   r13
  94956a:	push   r12
  94956c:	mov    r12,rdi
  94956f:	push   rbp
  949570:	mov    rbp,rcx
  949573:	push   rbx
  949574:	mov    rbx,rdx
  949577:	sub    rsp,0xa8
  94957e:	mov    QWORD PTR [rsp],rsi
  949582:	movzx  r15d,BYTE PTR [rip+0x2d5b58]        # c1f0e2 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  94958a:	test   r15b,r15b
  94958d:	je     949a34 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4d4>
  949593:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  94959f:	cmp    QWORD PTR [rax-0x1d8],0x0
  9495a7:	mov    r14,rax
  9495aa:	lea    rdx,[rax-0x360]
  9495b1:	je     9495dc <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x7c>
  9495b3:	movzx  r13d,BYTE PTR [rax-0x1e0]
  9495bb:	cmp    r13b,0xf
  9495bf:	jbe    949a9a <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x53a>
  9495c5:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  9495d1:	add    r13d,0x1
  9495d5:	mov    BYTE PTR [rax-0x1e0],r13b
  9495dc:	mov    rdx,QWORD PTR [rsp]
  9495e0:	mov    r8,QWORD PTR [r12+0x8]
  9495e5:	mov    rcx,QWORD PTR [r12]
  9495e9:	mov    rax,QWORD PTR [rdx+0x48]
  9495ed:	mov    rdi,QWORD PTR [rdx+0x38]
  9495f1:	mov    rdx,r8
  9495f4:	sub    rdx,rcx
  9495f7:	imul   rdi,rax
  9495fb:	sar    rdx,0x2
  9495ff:	lea    rsi,[rdi*4+0x0]
  949607:	cmp    rdx,rsi
  94960a:	jb     949a82 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x522>
  949610:	cmp    rsi,rdx
  949613:	jb     949a1a <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4ba>
  949619:	mov    rdx,QWORD PTR [rsp]
  94961d:	mov    r11,QWORD PTR [rdx+0x40]
  949621:	test   r11,r11
  949624:	je     9498d4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x374>
  94962a:	mov    r15,QWORD PTR [r12]
  94962e:	mov    rsi,QWORD PTR [rdx]
  949631:	xor    r12d,r12d
  949634:	lea    rdx,[rax*8+0x0]
  94963c:	lea    r13,[rax*4+0x0]
  949644:	shl    rax,0x4
  949648:	mov    QWORD PTR [rsp+0x68],r11
  94964d:	mov    rbx,QWORD PTR [rbx]
  949650:	lea    rcx,[rsi+rdx*1]
  949654:	add    rax,r15
  949657:	lea    r14,[r15+rdx*1]
  94965b:	mov    QWORD PTR [rsp+0x70],rsi
  949660:	mov    QWORD PTR [rsp+0x58],rax
  949665:	lea    rdi,[rcx+rdx*1]
  949669:	lea    rax,[r14+rdx*1]
  94966d:	mov    QWORD PTR [rsp+0x50],rdi
  949672:	mov    QWORD PTR [rsp+0x60],rax
  949677:	mov    QWORD PTR [rsp],0x0
  94967f:	mov    QWORD PTR [rsp+0x78],rbp
  949684:	mov    rbp,r13
  949687:	mov    r13,rcx
  94968a:	nop    WORD PTR [rax+rax*1+0x0]
  949690:	mov    rax,QWORD PTR [rsp+0x70]
  949695:	lea    rdi,[r15+r12*1]
  949699:	movss  xmm0,DWORD PTR [rax+r12*1]
  94969f:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  9496a4:	movaps xmm8,xmm0
  9496a8:	movss  xmm0,DWORD PTR [rax+rbp*1]
  9496ad:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  9496b2:	movaps xmm9,xmm0
  9496b6:	movss  xmm0,DWORD PTR [r13+r12*1+0x0]
  9496bd:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  9496c2:	movaps xmm10,xmm0
  9496c6:	movss  xmm0,DWORD PTR [r13+rbp*1+0x0]
  9496cd:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  9496d2:	mov    rax,QWORD PTR [rsp+0x50]
  9496d7:	movaps xmm11,xmm0
  9496db:	movss  xmm0,DWORD PTR [rax+r12*1]
  9496e1:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  9496e6:	movaps xmm12,xmm0
  9496ea:	movss  xmm0,DWORD PTR [rax+rbp*1]
  9496ef:	lea    rax,[r15+rbp*1]
  9496f3:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  9496f8:	movaps xmm13,xmm0
  9496fc:	movss  xmm0,DWORD PTR [r15+r12*1]
  949702:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949707:	mov    QWORD PTR [rsp+0x30],rax
  94970c:	lea    rax,[r14+r12*1]
  949710:	movaps xmm6,xmm0
  949713:	movss  xmm0,DWORD PTR [r15+rbp*1]
  949719:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94971e:	mov    QWORD PTR [rsp+0x28],rax
  949723:	lea    rax,[r14+rbp*1]
  949727:	movaps xmm5,xmm0
  94972a:	movss  xmm0,DWORD PTR [r14+r12*1]
  949730:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949735:	mov    QWORD PTR [rsp+0x10],rax
  94973a:	movaps xmm4,xmm0
  94973d:	movss  xmm0,DWORD PTR [r14+rbp*1]
  949743:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  949748:	mov    rcx,QWORD PTR [rsp+0x60]
  94974d:	mov    rax,QWORD PTR [rsp+0x58]
  949752:	movaps xmm3,xmm0
  949755:	movss  xmm0,DWORD PTR [rcx+r12*1]
  94975b:	lea    rdx,[rax+r12*1]
  94975f:	add    rax,rbp
  949762:	mov    QWORD PTR [rsp+0x18],rdx
  949767:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94976c:	mov    QWORD PTR [rsp+0x20],rax
  949771:	movaps xmm2,xmm0
  949774:	movss  xmm0,DWORD PTR [rcx+rbp*1]
  949779:	call   920600 <float cpl::simd::load<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  94977e:	mov    rsi,QWORD PTR [rsp+0x78]
  949783:	movaps xmm1,xmm0
  949786:	test   rsi,rsi
  949789:	je     949830 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x2d0>
  94978f:	xor    eax,eax
  949791:	nop    DWORD PTR [rax+0x0]
  949798:	movaps xmm7,xmm6
  94979b:	movaps xmm14,xmm5
  94979f:	movss  xmm0,DWORD PTR [rbx+rax*4]
  9497a4:	add    rax,0x1
  9497a8:	mulss  xmm14,xmm9
  9497ad:	call   9205f0 <float cpl::simd::broadcast<float>(cpl::simd::scalar_of<float, 4ul>::type const*) [clone .isra.0]>
  9497b2:	mulss  xmm7,xmm8
  9497b7:	mulss  xmm6,xmm9
  9497bc:	mulss  xmm5,xmm8
  9497c1:	subss  xmm7,xmm14
  9497c6:	movaps xmm14,xmm11
  9497ca:	mulss  xmm14,xmm3
  9497cf:	addss  xmm5,xmm6
  9497d3:	mulss  xmm3,xmm10
  9497d8:	movaps xmm6,xmm7
  9497db:	movaps xmm7,xmm4
  9497de:	mulss  xmm7,xmm10
  9497e3:	addss  xmm6,xmm0
  9497e7:	mulss  xmm4,xmm11
  9497ec:	subss  xmm7,xmm14
  9497f1:	movaps xmm14,xmm1
  9497f5:	addss  xmm3,xmm4
  9497f9:	mulss  xmm14,xmm13
  9497fe:	mulss  xmm1,xmm12
  949803:	movaps xmm4,xmm7
  949806:	movaps xmm7,xmm2
  949809:	mulss  xmm7,xmm12
  94980e:	addss  xmm4,xmm0
  949812:	mulss  xmm2,xmm13
  949817:	subss  xmm7,xmm14
  94981c:	addss  xmm1,xmm2
  949820:	movaps xmm2,xmm0
  949823:	addss  xmm2,xmm7
  949827:	cmp    rsi,rax
  94982a:	jne    949798 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x238>
  949830:	movaps xmm0,xmm6
  949833:	movss  DWORD PTR [rsp+0x3c],xmm1
  949839:	add    r12,0x4
  94983d:	add    rbp,0x4
  949841:	movss  DWORD PTR [rsp+0x40],xmm2
  949847:	movss  DWORD PTR [rsp+0x48],xmm3
  94984d:	movss  DWORD PTR [rsp+0x4c],xmm4
  949853:	movss  DWORD PTR [rsp+0x44],xmm5
  949859:	call   85fb30 <cpl::simd::store(float*, float)>
  94985e:	movss  xmm5,DWORD PTR [rsp+0x44]
  949864:	mov    rdi,QWORD PTR [rsp+0x30]
  949869:	movaps xmm0,xmm5
  94986c:	call   85fb30 <cpl::simd::store(float*, float)>
  949871:	movss  xmm4,DWORD PTR [rsp+0x4c]
  949877:	mov    rdi,QWORD PTR [rsp+0x28]
  94987c:	movaps xmm0,xmm4
  94987f:	call   85fb30 <cpl::simd::store(float*, float)>
  949884:	movss  xmm3,DWORD PTR [rsp+0x48]
  94988a:	mov    rdi,QWORD PTR [rsp+0x10]
  94988f:	movaps xmm0,xmm3
  949892:	call   85fb30 <cpl::simd::store(float*, float)>
  949897:	movss  xmm2,DWORD PTR [rsp+0x40]
  94989d:	mov    rdi,QWORD PTR [rsp+0x18]
  9498a2:	movaps xmm0,xmm2
  9498a5:	call   85fb30 <cpl::simd::store(float*, float)>
  9498aa:	movss  xmm1,DWORD PTR [rsp+0x3c]
  9498b0:	mov    rdi,QWORD PTR [rsp+0x20]
  9498b5:	movaps xmm0,xmm1
  9498b8:	call   85fb30 <cpl::simd::store(float*, float)>
  9498bd:	add    QWORD PTR [rsp],0x1
  9498c2:	mov    rdx,QWORD PTR [rsp+0x68]
  9498c7:	mov    rax,QWORD PTR [rsp]
  9498cb:	cmp    rax,rdx
  9498ce:	jne    949690 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x130>
  9498d4:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  9498e0:	lea    rbp,[rax-0x360]
  9498e7:	mov    r12,rax
  9498ea:	mov    r13,QWORD PTR [rbp+0x188]
  9498f1:	test   r13,r13
  9498f4:	je     949914 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3b4>
  9498f6:	movzx  eax,BYTE PTR [rbp+0x180]
  9498fd:	lea    ebx,[rax-0x1]
  949900:	mov    BYTE PTR [rbp+0x180],bl
  949906:	cmp    bl,0xf
  949909:	jbe    949926 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3c6>
  94990b:	add    WORD PTR [r13+0xfe2],0x1
  949914:	add    rsp,0xa8
  94991b:	pop    rbx
  94991c:	pop    rbp
  94991d:	pop    r12
  94991f:	pop    r13
  949921:	pop    r14
  949923:	pop    r15
  949925:	ret
  949926:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  94992b:	movzx  r14d,bl
  94992f:	lea    rdx,[r14+r14*2]
  949933:	mov    rcx,rax
  949936:	movq   xmm0,rax
  94993b:	shl    rdx,0x3
  94993f:	movdqu xmm1,XMMWORD PTR [r12+rdx*1-0x360]
  949949:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  94994e:	movq   xmm3,rcx
  949953:	movzx  r12d,BYTE PTR [rbp+0x181]
  94995b:	punpcklqdq xmm0,xmm3
  94995f:	psubq  xmm0,xmm1
  949963:	cmp    r12b,bl
  949966:	jae    949986 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x426>
  949968:	lea    eax,[r14-0x1]
  94996c:	movdqa xmm1,xmm0
  949970:	cdqe
  949972:	lea    rax,[rax+rax*2]
  949976:	movq   xmm2,QWORD PTR [rbp+rax*8+0x8]
  94997c:	paddq  xmm1,xmm2
  949980:	movq   QWORD PTR [rbp+rax*8+0x8],xmm1
  949986:	movzx  ebp,WORD PTR [r13+0xfe0]
  94998e:	cmp    bp,0x7f
  949992:	je     949af2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x592>
  949998:	movaps XMMWORD PTR [rsp],xmm0
  94999c:	sub    ebx,r12d
  94999f:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  9499ab:	lea    rdx,[r14+r14*2]
  9499af:	movdqa xmm0,XMMWORD PTR [rsp]
  9499b4:	mov    DWORD PTR [rsp+0x98],0x0
  9499bf:	movups XMMWORD PTR [rsp+0x88],xmm0
  9499c7:	lea    rdx,[rax+rdx*8-0x360]
  9499cf:	mov    rax,QWORD PTR [rdx]
  9499d2:	mov    QWORD PTR [rsp+0x80],rax
  9499da:	movzx  eax,BYTE PTR [rdx+0x10]
  9499de:	lea    edx,[rbp+0x1]
  9499e1:	shl    rbp,0x5
  9499e5:	movdqa xmm5,XMMWORD PTR [rsp+0x80]
  9499ee:	mov    WORD PTR [r13+0xfe0],dx
  9499f6:	mov    ah,bl
  9499f8:	mov    WORD PTR [rsp+0x9c],ax
  949a00:	lea    rax,[r13+rbp*1+0x0]
  949a05:	movdqu xmm3,XMMWORD PTR [rsp+0x8e]
  949a0e:	movups XMMWORD PTR [rax],xmm5
  949a11:	movups XMMWORD PTR [rax+0xe],xmm3
  949a15:	jmp    949914 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3b4>
  949a1a:	shl    rdi,0x4
  949a1e:	add    rcx,rdi
  949a21:	cmp    r8,rcx
  949a24:	je     949619 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb9>
  949a2a:	mov    QWORD PTR [r12+0x8],rcx
  949a2f:	jmp    949619 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb9>
  949a34:	mov    eax,0x1
  949a39:	lock xadd DWORD PTR [rip+0x2ca11f],eax        # c13b60 <cpl::Profiling::registerRegion(char const*)::counter>
  949a41:	mov    edx,0x1
  949a46:	mov    r15d,0x1
  949a4c:	add    eax,0x2
  949a4f:	cmp    eax,0xfe
  949a54:	ja     949a6f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x50f>
  949a56:	lea    rdx,[rip+0x2c85c3]        # c12020 <cpl::Profiling::regions>
  949a5d:	mov    ecx,eax
  949a5f:	lea    rdi,[rip+0x9e6e0]        # 9e8146 <_IO_stdin_used+0x18146>
  949a66:	mov    r15d,eax
  949a69:	mov    QWORD PTR [rdx+rcx*8],rdi
  949a6d:	mov    edx,eax
  949a6f:	xor    eax,eax
  949a71:	lock cmpxchg BYTE PTR [rip+0x2d5669],dl        # c1f0e2 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  949a79:	cmovne r15d,eax
  949a7d:	jmp    949593 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x33>
  949a82:	sub    rsi,rdx
  949a85:	mov    rdi,r12
  949a88:	call   91d820 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  949a8d:	mov    rax,QWORD PTR [rsp]
  949a91:	mov    rax,QWORD PTR [rax+0x48]
  949a95:	jmp    949619 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb9>
  949a9a:	mov    QWORD PTR [rsp+0x10],rdx
  949a9f:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  949aa4:	mov    rdx,QWORD PTR [rsp+0x10]
  949aa9:	mov    QWORD PTR [rsp+0x88],0x0
  949ab5:	mov    QWORD PTR [rsp+0x80],rax
  949abd:	lea    rax,[r13+r13*2+0x0]
  949ac2:	movdqa xmm5,XMMWORD PTR [rsp+0x80]
  949acb:	shl    rax,0x3
  949acf:	mov    BYTE PTR [rsp+0x90],r15b
  949ad7:	movups XMMWORD PTR [r14+rax*1-0x360],xmm5
  949ae0:	mov    BYTE PTR [rax+rdx*1+0x10],r15b
  949ae5:	movzx  r13d,BYTE PTR [rdx+0x180]
  949aed:	jmp    9495c5 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x65>
  949af2:	movzx  eax,WORD PTR [r13+0xfe2]
  949afa:	cmp    ax,0xffff
  949afe:	je     949914 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3b4>
  949b04:	add    eax,0x1
  949b07:	mov    WORD PTR [r13+0xfe2],ax
  949b0f:	jmp    949914 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 1ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3b4>
  949b14:	endbr64
  949b18:	mov    rbx,rax
  949b1b:	call   920350 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  949b20:	mov    rdi,rbx
  949b23:	call   e1b80 <_Unwind_Resume@plt>

Disassembly of section .fini:
