; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000931010 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long):
  931010:	endbr64
  931014:	lea    r10,[rsp+0x8]
  931019:	and    rsp,0xffffffffffffffe0
  93101d:	push   QWORD PTR [r10-0x8]
  931021:	push   rbp
  931022:	mov    rbp,rsp
  931025:	push   r15
  931027:	push   r14
  931029:	push   r13
  93102b:	push   r12
  93102d:	push   r10
  93102f:	push   rbx
  931030:	sub    rsp,0x740
  931037:	mov    QWORD PTR [rbp-0x698],rdi
  93103e:	mov    QWORD PTR [rbp-0x6a8],rsi
  931045:	mov    QWORD PTR [rbp-0x6b0],rdx
  93104c:	mov    QWORD PTR [rbp-0x650],rcx
  931053:	mov    rax,QWORD PTR fs:0x28
  93105c:	mov    QWORD PTR [rbp-0x38],rax
  931060:	xor    eax,eax
  931062:	movzx  r13d,BYTE PTR [rip+0x2ee09e]        # c1f108 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  93106a:	test   r13b,r13b
  93106d:	je     931ce7 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xcd7>
  931073:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  93107f:	cmp    QWORD PTR [rax-0x1d8],0x0
  931087:	mov    rbx,rax
  93108a:	lea    r14,[rax-0x360]
  931091:	je     9310bc <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xac>
  931093:	movzx  r12d,BYTE PTR [r14+0x180]
  93109b:	cmp    r12b,0xf
  93109f:	jbe    931d58 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xd48>
  9310a5:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  9310b1:	add    r12d,0x1
  9310b5:	mov    BYTE PTR [rax-0x1e0],r12b
  9310bc:	mov    rdx,QWORD PTR [rbp-0x6a8]
  9310c3:	mov    rcx,QWORD PTR [rbp-0x698]
  9310ca:	mov    rax,QWORD PTR [rdx+0x48]
  9310ce:	mov    rdx,QWORD PTR [rdx+0x38]
  9310d2:	mov    r9,QWORD PTR [rcx+0x8]
  9310d6:	mov    rdi,QWORD PTR [rcx]
  9310d9:	mov    r8,rdx
  9310dc:	imul   r8,rax
  9310e0:	mov    rcx,r9
  9310e3:	sub    rcx,rdi
  9310e6:	sar    rcx,0x2
  9310ea:	lea    rsi,[r8*4+0x0]
  9310f2:	cmp    rcx,rsi
  9310f5:	jb     931d35 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xd25>
  9310fb:	cmp    rsi,rcx
  9310fe:	jb     931cc7 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xcb7>
  931104:	mov    rdi,QWORD PTR [rbp-0x6a8]
  93110b:	lea    rcx,[rax+rax*1]
  93110f:	imul   rdx,rcx
  931113:	cmp    QWORD PTR [rdi+0x40],0x0
  931118:	je     931b4e <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb3e>
  93111e:	lea    rdi,[rcx+rax*1]
  931122:	lea    r10,[rax*4+0x0]
  93112a:	mov    QWORD PTR [rbp-0x680],0x0
  931135:	shl    rdi,0x2
  931139:	lea    rsi,[r10+rdx*1]
  93113d:	mov    QWORD PTR [rbp-0x678],r10
  931144:	lea    r13,[rbp-0xd0]
  93114b:	mov    QWORD PTR [rbp-0x670],rdi
  931152:	mov    rdi,rax
  931155:	add    rax,r10
  931158:	lea    r9,[rbp-0x1d0]
  93115f:	shl    rax,0x2
  931163:	shl    rdi,0x4
  931167:	mov    QWORD PTR [rbp-0x6a0],0x0
  931172:	lea    r8,[rbp-0x190]
  931179:	mov    QWORD PTR [rbp-0x660],rax
  931180:	lea    rax,[rdx*4+0x0]
  931188:	lea    r12,[rbp-0x90]
  93118f:	mov    QWORD PTR [rbp-0x690],rax
  931196:	mov    rax,rsi
  931199:	lea    rdx,[rbp-0x110]
  9311a0:	sub    rax,rcx
  9311a3:	mov    QWORD PTR [rbp-0x668],rdi
  9311aa:	lea    rcx,[rbp-0x150]
  9311b1:	shl    rax,0x2
  9311b5:	mov    QWORD PTR [rbp-0x758],r10
  9311bc:	mov    r14,rcx
  9311bf:	mov    QWORD PTR [rbp-0x6c0],rax
  9311c6:	lea    rax,[rsi*4+0x0]
  9311ce:	mov    QWORD PTR [rbp-0x6b8],rax
  9311d5:	lea    rax,[rbp-0x470]
  9311dc:	mov    QWORD PTR [rbp-0x6d8],rax
  9311e3:	lea    rax,[rbp-0x450]
  9311ea:	mov    QWORD PTR [rbp-0x6f0],rax
  9311f1:	lea    rax,[rbp-0x430]
  9311f8:	mov    QWORD PTR [rbp-0x708],rax
  9311ff:	lea    rax,[rbp-0x410]
  931206:	mov    QWORD PTR [rbp-0x720],rax
  93120d:	lea    rax,[rbp-0x3f0]
  931214:	mov    QWORD PTR [rbp-0x738],rax
  93121b:	lea    rax,[rbp-0x3d0]
  931222:	mov    QWORD PTR [rbp-0x750],rax
  931229:	lea    rax,[rbp-0x390]
  931230:	mov    QWORD PTR [rbp-0x6e8],rax
  931237:	lea    rax,[rbp-0x370]
  93123e:	mov    QWORD PTR [rbp-0x700],rax
  931245:	lea    rax,[rbp-0x350]
  93124c:	mov    QWORD PTR [rbp-0x718],rax
  931253:	lea    rax,[rbp-0x330]
  93125a:	mov    QWORD PTR [rbp-0x730],rax
  931261:	lea    rax,[rbp-0x310]
  931268:	mov    QWORD PTR [rbp-0x748],rax
  93126f:	lea    rax,[rbp-0x2f0]
  931276:	mov    QWORD PTR [rbp-0x6d0],rax
  93127d:	lea    rax,[rbp-0x2d0]
  931284:	mov    QWORD PTR [rbp-0x6e0],rax
  93128b:	lea    rax,[rbp-0x2b0]
  931292:	mov    QWORD PTR [rbp-0x6f8],rax
  931299:	lea    rax,[rbp-0x290]
  9312a0:	mov    QWORD PTR [rbp-0x710],rax
  9312a7:	lea    rax,[rbp-0x270]
  9312ae:	mov    QWORD PTR [rbp-0x728],rax
  9312b5:	lea    rax,[rbp-0x250]
  9312bc:	mov    QWORD PTR [rbp-0x740],rax
  9312c3:	lea    rax,[rbp-0x230]
  9312ca:	mov    QWORD PTR [rbp-0x6c8],rax
  9312d1:	lea    rax,[rbp-0x3b0]
  9312d8:	mov    QWORD PTR [rbp-0x640],rax
  9312df:	mov    rax,r13
  9312e2:	nop    WORD PTR [rax+rax*1+0x0]
  9312e8:	mov    QWORD PTR [rbp-0x630],rax
  9312ef:	mov    rax,QWORD PTR [rbp-0x6a8]
  9312f6:	mov    r13,QWORD PTR [rbp-0x680]
  9312fd:	mov    rdi,QWORD PTR [rbp-0x6d8]
  931304:	mov    QWORD PTR [rbp-0x638],r8
  93130b:	mov    rbx,QWORD PTR [rax]
  93130e:	mov    QWORD PTR [rbp-0x628],r9
  931315:	mov    QWORD PTR [rbp-0x620],rdx
  93131c:	lea    rsi,[rbx+r13*1]
  931320:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  931325:	mov    r15,QWORD PTR [rbp-0x678]
  93132c:	mov    rdi,QWORD PTR [rbp-0x6f0]
  931333:	lea    rsi,[rbx+r15*1]
  931337:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  93133c:	mov    rax,QWORD PTR [rbp-0x758]
  931343:	mov    rdi,QWORD PTR [rbp-0x708]
  93134a:	add    rax,r15
  93134d:	lea    rsi,[rbx+rax*1]
  931351:	mov    QWORD PTR [rbp-0x688],rax
  931358:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  93135d:	mov    rdx,QWORD PTR [rbp-0x670]
  931364:	mov    rdi,QWORD PTR [rbp-0x720]
  93136b:	lea    rsi,[rbx+rdx*1]
  93136f:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  931374:	mov    rcx,QWORD PTR [rbp-0x668]
  93137b:	mov    rdi,QWORD PTR [rbp-0x738]
  931382:	lea    rsi,[rbx+rcx*1]
  931386:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  93138b:	mov    rdi,QWORD PTR [rbp-0x660]
  931392:	lea    rsi,[rbx+rdi*1]
  931396:	mov    rdi,QWORD PTR [rbp-0x750]
  93139d:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9313a2:	mov    rbx,QWORD PTR [rbp-0x6b0]
  9313a9:	mov    rdi,QWORD PTR [rbp-0x6e8]
  9313b0:	mov    rsi,QWORD PTR [rbx]
  9313b3:	mov    QWORD PTR [rbp-0x1e0],rsi
  9313ba:	mov    rsi,QWORD PTR [rbp-0x698]
  9313c1:	mov    rbx,QWORD PTR [rsi]
  9313c4:	lea    rsi,[rbx+r13*1]
  9313c8:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9313cd:	mov    rdi,QWORD PTR [rbp-0x700]
  9313d4:	lea    rsi,[rbx+r15*1]
  9313d8:	movdqa xmm4,XMMWORD PTR [rbp-0x390]
  9313e0:	movdqa xmm3,XMMWORD PTR [rbp-0x380]
  9313e8:	movaps XMMWORD PTR [rbp-0x1d0],xmm4
  9313ef:	movaps XMMWORD PTR [rbp-0x1c0],xmm3
  9313f6:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9313fb:	mov    rax,QWORD PTR [rbp-0x688]
  931402:	movdqa xmm7,XMMWORD PTR [rbp-0x370]
  93140a:	movdqa xmm2,XMMWORD PTR [rbp-0x360]
  931412:	mov    rdi,QWORD PTR [rbp-0x718]
  931419:	lea    rsi,[rbx+rax*1]
  93141d:	movaps XMMWORD PTR [rbp-0x190],xmm7
  931424:	movaps XMMWORD PTR [rbp-0x180],xmm2
  93142b:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  931430:	mov    rdx,QWORD PTR [rbp-0x670]
  931437:	movdqa xmm4,XMMWORD PTR [rbp-0x350]
  93143f:	movdqa xmm3,XMMWORD PTR [rbp-0x340]
  931447:	mov    rdi,QWORD PTR [rbp-0x730]
  93144e:	lea    rsi,[rbx+rdx*1]
  931452:	movaps XMMWORD PTR [rbp-0x150],xmm4
  931459:	movaps XMMWORD PTR [rbp-0x140],xmm3
  931460:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  931465:	mov    rcx,QWORD PTR [rbp-0x668]
  93146c:	movdqa xmm7,XMMWORD PTR [rbp-0x330]
  931474:	movdqa xmm2,XMMWORD PTR [rbp-0x320]
  93147c:	mov    rdi,QWORD PTR [rbp-0x748]
  931483:	lea    rsi,[rbx+rcx*1]
  931487:	movaps XMMWORD PTR [rbp-0x110],xmm7
  93148e:	movaps XMMWORD PTR [rbp-0x100],xmm2
  931495:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  93149a:	mov    rdi,QWORD PTR [rbp-0x660]
  9314a1:	movdqa xmm4,XMMWORD PTR [rbp-0x310]
  9314a9:	movdqa xmm3,XMMWORD PTR [rbp-0x300]
  9314b1:	lea    rsi,[rbx+rdi*1]
  9314b5:	mov    rdi,QWORD PTR [rbp-0x6d0]
  9314bc:	movaps XMMWORD PTR [rbp-0xd0],xmm4
  9314c3:	movaps XMMWORD PTR [rbp-0xc0],xmm3
  9314ca:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9314cf:	mov    rax,QWORD PTR [rbp-0x6b0]
  9314d6:	movdqa xmm7,XMMWORD PTR [rbp-0x2f0]
  9314de:	movdqa xmm2,XMMWORD PTR [rbp-0x2e0]
  9314e6:	mov    rdi,QWORD PTR [rbp-0x6e0]
  9314ed:	mov    rsi,QWORD PTR [rax+0x10]
  9314f1:	mov    rax,QWORD PTR [rbp-0x690]
  9314f8:	movaps XMMWORD PTR [rbp-0x90],xmm7
  9314ff:	movaps XMMWORD PTR [rbp-0x80],xmm2
  931503:	mov    QWORD PTR [rbp-0x1d8],rsi
  93150a:	lea    rsi,[rax+r13*1]
  93150e:	add    rsi,rbx
  931511:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  931516:	mov    rax,QWORD PTR [rbp-0x690]
  93151d:	movdqa xmm4,XMMWORD PTR [rbp-0x2d0]
  931525:	movdqa xmm3,XMMWORD PTR [rbp-0x2c0]
  93152d:	mov    rdi,QWORD PTR [rbp-0x6f8]
  931534:	lea    rsi,[rax+r15*1]
  931538:	movaps XMMWORD PTR [rbp-0x1b0],xmm4
  93153f:	add    rsi,rbx
  931542:	movaps XMMWORD PTR [rbp-0x1a0],xmm3
  931549:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  93154e:	mov    rax,QWORD PTR [rbp-0x6c0]
  931555:	movdqa xmm7,XMMWORD PTR [rbp-0x2b0]
  93155d:	movdqa xmm2,XMMWORD PTR [rbp-0x2a0]
  931565:	mov    rdi,QWORD PTR [rbp-0x710]
  93156c:	lea    rsi,[rax+r13*1]
  931570:	movaps XMMWORD PTR [rbp-0x170],xmm7
  931577:	add    rsi,rbx
  93157a:	movaps XMMWORD PTR [rbp-0x160],xmm2
  931581:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  931586:	mov    rax,QWORD PTR [rbp-0x6c0]
  93158d:	movdqa xmm3,XMMWORD PTR [rbp-0x280]
  931595:	movdqa xmm4,XMMWORD PTR [rbp-0x290]
  93159d:	mov    rdi,QWORD PTR [rbp-0x728]
  9315a4:	lea    rsi,[rax+r15*1]
  9315a8:	movaps XMMWORD PTR [rbp-0x120],xmm3
  9315af:	add    rsi,rbx
  9315b2:	movaps XMMWORD PTR [rbp-0x130],xmm4
  9315b9:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9315be:	mov    rax,QWORD PTR [rbp-0x6b8]
  9315c5:	movdqa xmm7,XMMWORD PTR [rbp-0x270]
  9315cd:	movdqa xmm2,XMMWORD PTR [rbp-0x260]
  9315d5:	mov    rdi,QWORD PTR [rbp-0x740]
  9315dc:	lea    rsi,[rax+r13*1]
  9315e0:	movaps XMMWORD PTR [rbp-0xf0],xmm7
  9315e7:	add    rsi,rbx
  9315ea:	movaps XMMWORD PTR [rbp-0xe0],xmm2
  9315f1:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9315f6:	mov    rax,QWORD PTR [rbp-0x6b8]
  9315fd:	movdqa xmm5,XMMWORD PTR [rbp-0x250]
  931605:	movdqa xmm4,XMMWORD PTR [rbp-0x240]
  93160d:	mov    rdi,QWORD PTR [rbp-0x6c8]
  931614:	lea    rsi,[rax+r15*1]
  931618:	movaps XMMWORD PTR [rbp-0xb0],xmm5
  93161f:	add    rsi,rbx
  931622:	movaps XMMWORD PTR [rbp-0xa0],xmm4
  931629:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  93162e:	movdqa xmm3,XMMWORD PTR [rbp-0x230]
  931636:	xor    esi,esi
  931638:	movdqa xmm7,XMMWORD PTR [rbp-0x220]
  931640:	cmp    QWORD PTR [rbp-0x650],0x0
  931648:	mov    rdx,QWORD PTR [rbp-0x620]
  93164f:	lea    rdi,[rbp-0x1e0]
  931656:	mov    r9,QWORD PTR [rbp-0x628]
  93165d:	mov    rax,QWORD PTR [rbp-0x630]
  931664:	movaps XMMWORD PTR [rbp-0x70],xmm3
  931668:	mov    r8,QWORD PTR [rbp-0x638]
  93166f:	movaps XMMWORD PTR [rbp-0x60],xmm7
  931673:	je     931966 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x956>
  931679:	mov    QWORD PTR [rbp-0x658],rdi
  931680:	mov    r15,r14
  931683:	mov    r14,r12
  931686:	mov    r12,rdx
  931689:	mov    rdx,rsi
  93168c:	nop    DWORD PTR [rax+0x0]
  931690:	mov    rcx,QWORD PTR [rbp-0x658]
  931697:	mov    rdi,r15
  93169a:	mov    QWORD PTR [rbp-0x648],rdx
  9316a1:	mov    rdx,r12
  9316a4:	mov    r15,r14
  9316a7:	xor    ebx,ebx
  9316a9:	mov    r14,rdi
  9316ac:	mov    r12,rcx
  9316af:	mov    rcx,rdx
  9316b2:	mov    r13,QWORD PTR [r12]
  9316b6:	mov    rdi,QWORD PTR [rbp-0x640]
  9316bd:	mov    QWORD PTR [rbp-0x638],r8
  9316c4:	add    r12,0x8
  9316c8:	mov    QWORD PTR [rbp-0x630],rax
  9316cf:	mov    rsi,r13
  9316d2:	mov    QWORD PTR [rbp-0x628],r9
  9316d9:	add    r13,0x4
  9316dd:	mov    QWORD PTR [rbp-0x620],rcx
  9316e4:	call   9243d0 <float __vector(8) cpl::simd::broadcast<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9316e9:	mov    r9,QWORD PTR [rbp-0x628]
  9316f0:	mov    r8,QWORD PTR [rbp-0x638]
  9316f7:	movaps xmm0,XMMWORD PTR [rbp-0x450]
  9316fe:	movaps xmm2,XMMWORD PTR [rbp-0x470]
  931705:	lea    rsi,[r9+rbx*1]
  931709:	movaps xmm3,XMMWORD PTR [rbp-0x460]
  931710:	mov    rcx,QWORD PTR [rbp-0x620]
  931717:	movaps xmm1,XMMWORD PTR [rsi+0x10]
  93171b:	movaps xmm6,XMMWORD PTR [rsi]
  93171e:	lea    rsi,[r8+rbx*1]
  931722:	movaps xmm5,XMMWORD PTR [rsi]
  931725:	movaps xmm4,XMMWORD PTR [rsi+0x10]
  931729:	lea    rsi,[r14+rbx*1]
  93172d:	mulps  xmm0,xmm6
  931730:	mulps  xmm2,xmm5
  931733:	mulps  xmm3,xmm4
  931736:	mulps  xmm6,XMMWORD PTR [rbp-0x470]
  93173d:	mulps  xmm5,XMMWORD PTR [rbp-0x450]
  931744:	mulps  xmm4,XMMWORD PTR [rbp-0x440]
  93174b:	addps  xmm0,xmm2
  93174e:	movaps xmm2,XMMWORD PTR [rbp-0x440]
  931755:	mulps  xmm2,xmm1
  931758:	mulps  xmm1,XMMWORD PTR [rbp-0x460]
  93175f:	addps  xmm2,xmm3
  931762:	subps  xmm1,xmm4
  931765:	movaps XMMWORD PTR [rbp-0x600],xmm2
  93176c:	movaps XMMWORD PTR [r8+rbx*1],xmm0
  931771:	movaps xmm0,xmm6
  931774:	movdqa xmm7,XMMWORD PTR [rbp-0x600]
  93177c:	subps  xmm0,xmm5
  93177f:	movaps XMMWORD PTR [r8+rbx*1+0x10],xmm7
  931785:	movaps xmm3,XMMWORD PTR [rbp-0x3b0]
  93178c:	movaps xmm2,XMMWORD PTR [rbp-0x3a0]
  931793:	addps  xmm0,xmm3
  931796:	addps  xmm1,xmm2
  931799:	movaps XMMWORD PTR [rbp-0x610],xmm0
  9317a0:	movaps XMMWORD PTR [rbp-0x600],xmm1
  9317a7:	movaps xmm1,XMMWORD PTR [rbp-0x410]
  9317ae:	movaps XMMWORD PTR [r9+rbx*1],xmm0
  9317b3:	movdqa xmm4,XMMWORD PTR [rbp-0x600]
  9317bb:	movaps xmm7,XMMWORD PTR [rsi]
  9317be:	movaps xmm0,XMMWORD PTR [rbp-0x430]
  9317c5:	movaps XMMWORD PTR [r9+rbx*1+0x10],xmm4
  9317cb:	movaps xmm5,XMMWORD PTR [rsi+0x10]
  9317cf:	lea    rsi,[rcx+rbx*1]
  9317d3:	mulps  xmm1,xmm7
  9317d6:	movaps xmm6,XMMWORD PTR [rsi]
  9317d9:	movaps xmm4,XMMWORD PTR [rsi+0x10]
  9317dd:	mulps  xmm0,xmm6
  9317e0:	addps  xmm0,xmm1
  9317e3:	movaps xmm1,XMMWORD PTR [rbp-0x420]
  9317ea:	movaps xmm8,XMMWORD PTR [rbp-0x400]
  9317f2:	mulps  xmm7,XMMWORD PTR [rbp-0x430]
  9317f9:	mov    rax,QWORD PTR [rbp-0x630]
  931800:	mulps  xmm8,xmm5
  931804:	mulps  xmm1,xmm4
  931807:	lea    rsi,[rax+rbx*1]
  93180b:	mulps  xmm5,XMMWORD PTR [rbp-0x420]
  931812:	mulps  xmm6,XMMWORD PTR [rbp-0x410]
  931819:	mulps  xmm4,XMMWORD PTR [rbp-0x400]
  931820:	addps  xmm1,xmm8
  931824:	movaps xmm8,XMMWORD PTR [rbp-0x3c0]
  93182c:	movaps XMMWORD PTR [rbp-0x600],xmm1
  931833:	movaps xmm1,xmm5
  931836:	movaps XMMWORD PTR [rcx+rbx*1],xmm0
  93183a:	movdqa xmm0,XMMWORD PTR [rbp-0x600]
  931842:	subps  xmm1,xmm4
  931845:	movaps XMMWORD PTR [rcx+rbx*1+0x10],xmm0
  93184a:	movaps xmm0,xmm7
  93184d:	subps  xmm0,xmm6
  931850:	addps  xmm1,xmm2
  931853:	addps  xmm0,xmm3
  931856:	movaps XMMWORD PTR [rbp-0x600],xmm1
  93185d:	movaps xmm1,XMMWORD PTR [rbp-0x3d0]
  931864:	movaps XMMWORD PTR [rbp-0x610],xmm0
  93186b:	movaps XMMWORD PTR [r14+rbx*1],xmm0
  931870:	movdqa xmm7,XMMWORD PTR [rbp-0x600]
  931878:	movaps xmm0,XMMWORD PTR [rbp-0x3f0]
  93187f:	movaps XMMWORD PTR [r14+rbx*1+0x10],xmm7
  931885:	movaps xmm7,XMMWORD PTR [rsi]
  931888:	movaps xmm5,XMMWORD PTR [rsi+0x10]
  93188c:	lea    rsi,[r15+rbx*1]
  931890:	movaps xmm6,XMMWORD PTR [rsi]
  931893:	movaps xmm4,XMMWORD PTR [rsi+0x10]
  931897:	mulps  xmm1,xmm7
  93189a:	mulps  xmm0,xmm6
  93189d:	mulps  xmm8,xmm5
  9318a1:	mulps  xmm7,XMMWORD PTR [rbp-0x3f0]
  9318a8:	mulps  xmm6,XMMWORD PTR [rbp-0x3d0]
  9318af:	mulps  xmm5,XMMWORD PTR [rbp-0x3e0]
  9318b6:	addps  xmm0,xmm1
  9318b9:	movaps xmm1,XMMWORD PTR [rbp-0x3e0]
  9318c0:	mulps  xmm1,xmm4
  9318c3:	addps  xmm1,xmm8
  9318c7:	movaps XMMWORD PTR [rbp-0x600],xmm1
  9318ce:	movaps xmm1,xmm5
  9318d1:	movaps XMMWORD PTR [r15+rbx*1],xmm0
  9318d6:	movdqa xmm0,XMMWORD PTR [rbp-0x600]
  9318de:	movaps XMMWORD PTR [r15+rbx*1+0x10],xmm0
  9318e4:	movaps xmm0,xmm7
  9318e7:	subps  xmm0,xmm6
  9318ea:	addps  xmm0,xmm3
  9318ed:	movaps XMMWORD PTR [rbp-0x610],xmm0
  9318f4:	movaps xmm3,XMMWORD PTR [rbp-0x3c0]
  9318fb:	mulps  xmm3,xmm4
  9318fe:	subps  xmm1,xmm3
  931901:	addps  xmm1,xmm2
  931904:	movaps XMMWORD PTR [rbp-0x600],xmm1
  93190b:	movaps XMMWORD PTR [rax+rbx*1],xmm0
  93190f:	movdqa xmm3,XMMWORD PTR [rbp-0x600]
  931917:	mov    QWORD PTR [r12-0x8],r13
  93191c:	movaps XMMWORD PTR [rax+rbx*1+0x10],xmm3
  931921:	add    rbx,0x20
  931925:	cmp    rbx,0x40
  931929:	jne    9316b2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x6a2>
  93192f:	mov    rdx,QWORD PTR [rbp-0x648]
  931936:	mov    rdi,r14
  931939:	mov    r12,rcx
  93193c:	mov    r14,r15
  93193f:	mov    r15,rdi
  931942:	add    rdx,0x1
  931946:	cmp    QWORD PTR [rbp-0x650],rdx
  93194d:	jne    931690 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x680>
  931953:	mov    rdx,rcx
  931956:	mov    rcx,QWORD PTR [rbp-0x698]
  93195d:	mov    r12,r14
  931960:	mov    r14,rdi
  931963:	mov    rbx,QWORD PTR [rcx]
  931966:	xor    r15d,r15d
  931969:	mov    rdi,rbx
  93196c:	mov    rbx,r14
  93196f:	mov    r14,r12
  931972:	mov    r12,r15
  931975:	mov    r15,QWORD PTR [rbp-0x698]
  93197c:	xor    r13d,r13d
  93197f:	movdqa xmm4,XMMWORD PTR [r9+r12*1]
  931985:	sub    rsp,0x20
  931989:	mov    QWORD PTR [rbp-0x648],rax
  931990:	mov    rax,QWORD PTR [rbp-0x680]
  931997:	mov    QWORD PTR [rbp-0x620],rdx
  93199e:	mov    QWORD PTR [rbp-0x628],r8
  9319a5:	movaps XMMWORD PTR [rsp],xmm4
  9319a9:	movdqa xmm4,XMMWORD PTR [r9+r12*1+0x10]
  9319b0:	lea    rsi,[rax+r13*1]
  9319b4:	mov    QWORD PTR [rbp-0x638],r9
  9319bb:	add    rdi,rsi
  9319be:	movaps XMMWORD PTR [rsp+0x10],xmm4
  9319c3:	call   85fb10 <cpl::simd::store(float*, float __vector(8))>
  9319c8:	mov    r8,QWORD PTR [rbp-0x628]
  9319cf:	mov    rax,QWORD PTR [rbp-0x678]
  9319d6:	movdqa xmm4,XMMWORD PTR [r8+r12*1]
  9319dc:	lea    rdi,[rax+r13*1]
  9319e0:	add    rdi,QWORD PTR [r15]
  9319e3:	movaps XMMWORD PTR [rsp],xmm4
  9319e7:	movdqa xmm4,XMMWORD PTR [r8+r12*1+0x10]
  9319ee:	mov    QWORD PTR [rbp-0x630],r8
  9319f5:	movaps XMMWORD PTR [rsp+0x10],xmm4
  9319fa:	call   85fb10 <cpl::simd::store(float*, float __vector(8))>
  9319ff:	movdqa xmm4,XMMWORD PTR [rbx+r12*1]
  931a05:	mov    rax,QWORD PTR [rbp-0x688]
  931a0c:	lea    rdi,[rax+r13*1]
  931a10:	add    rdi,QWORD PTR [r15]
  931a13:	movaps XMMWORD PTR [rsp],xmm4
  931a17:	movdqa xmm3,XMMWORD PTR [rbx+r12*1+0x10]
  931a1e:	movaps XMMWORD PTR [rsp+0x10],xmm3
  931a23:	call   85fb10 <cpl::simd::store(float*, float __vector(8))>
  931a28:	mov    rdx,QWORD PTR [rbp-0x620]
  931a2f:	mov    rax,QWORD PTR [rbp-0x670]
  931a36:	movdqa xmm7,XMMWORD PTR [rdx+r12*1]
  931a3c:	lea    rdi,[rax+r13*1]
  931a40:	add    rdi,QWORD PTR [r15]
  931a43:	movaps XMMWORD PTR [rsp],xmm7
  931a47:	movdqa xmm2,XMMWORD PTR [rdx+r12*1+0x10]
  931a4e:	mov    QWORD PTR [rbp-0x628],rdx
  931a55:	movaps XMMWORD PTR [rsp+0x10],xmm2
  931a5a:	call   85fb10 <cpl::simd::store(float*, float __vector(8))>
  931a5f:	mov    rax,QWORD PTR [rbp-0x668]
  931a66:	lea    rdi,[rax+r13*1]
  931a6a:	mov    rax,QWORD PTR [rbp-0x648]
  931a71:	add    rdi,QWORD PTR [r15]
  931a74:	movdqa xmm4,XMMWORD PTR [rax+r12*1]
  931a7a:	movaps XMMWORD PTR [rsp],xmm4
  931a7e:	movdqa xmm3,XMMWORD PTR [rax+r12*1+0x10]
  931a85:	mov    QWORD PTR [rbp-0x620],rax
  931a8c:	movaps XMMWORD PTR [rsp+0x10],xmm3
  931a91:	call   85fb10 <cpl::simd::store(float*, float __vector(8))>
  931a96:	mov    rax,QWORD PTR [rbp-0x660]
  931a9d:	movdqa xmm7,XMMWORD PTR [r14+r12*1]
  931aa3:	lea    rdi,[rax+r13*1]
  931aa7:	add    rdi,QWORD PTR [r15]
  931aaa:	movaps XMMWORD PTR [rsp],xmm7
  931aae:	movdqa xmm2,XMMWORD PTR [r14+r12*1+0x10]
  931ab5:	add    r12,0x20
  931ab9:	movaps XMMWORD PTR [rsp+0x10],xmm2
  931abe:	call   85fb10 <cpl::simd::store(float*, float __vector(8))>
  931ac3:	add    rsp,0x20
  931ac7:	cmp    r12,0x40
  931acb:	mov    rax,QWORD PTR [rbp-0x620]
  931ad2:	mov    rdx,QWORD PTR [rbp-0x628]
  931ad9:	mov    r8,QWORD PTR [rbp-0x630]
  931ae0:	mov    r9,QWORD PTR [rbp-0x638]
  931ae7:	je     931b00 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xaf0>
  931ae9:	mov    rcx,QWORD PTR [rbp-0x690]
  931af0:	mov    rdi,QWORD PTR [r15]
  931af3:	add    r13,rcx
  931af6:	jmp    93197f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x96f>
  931afb:	nop    DWORD PTR [rax+rax*1+0x0]
  931b00:	mov    rdi,QWORD PTR [rbp-0x6a8]
  931b07:	add    QWORD PTR [rbp-0x6a0],0x8
  931b0f:	mov    r12,r14
  931b12:	mov    r14,rbx
  931b15:	add    QWORD PTR [rbp-0x680],0x20
  931b1d:	mov    rcx,QWORD PTR [rbp-0x6a0]
  931b24:	add    QWORD PTR [rbp-0x678],0x20
  931b2c:	add    QWORD PTR [rbp-0x670],0x20
  931b34:	add    QWORD PTR [rbp-0x668],0x20
  931b3c:	add    QWORD PTR [rbp-0x660],0x20
  931b44:	cmp    rcx,QWORD PTR [rdi+0x40]
  931b48:	jb     9312e8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x2d8>
  931b4e:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  931b5a:	lea    r12,[rax-0x360]
  931b61:	mov    r13,rax
  931b64:	mov    r14,QWORD PTR [r12+0x188]
  931b6c:	test   r14,r14
  931b6f:	je     931b93 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb83>
  931b71:	movzx  eax,BYTE PTR [r12+0x180]
  931b7a:	lea    ebx,[rax-0x1]
  931b7d:	mov    BYTE PTR [r12+0x180],bl
  931b85:	cmp    bl,0xf
  931b88:	jbe    931bbb <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xbab>
  931b8a:	add    WORD PTR [r14+0xfe2],0x1
  931b93:	mov    rax,QWORD PTR [rbp-0x38]
  931b97:	sub    rax,QWORD PTR fs:0x28
  931ba0:	jne    931dd6 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xdc6>
  931ba6:	lea    rsp,[rbp-0x30]
  931baa:	pop    rbx
  931bab:	pop    r10
  931bad:	pop    r12
  931baf:	pop    r13
  931bb1:	pop    r14
  931bb3:	pop    r15
  931bb5:	pop    rbp
  931bb6:	lea    rsp,[r10-0x8]
  931bba:	ret
  931bbb:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  931bc0:	movzx  r15d,bl
  931bc4:	movsxd rcx,r15d
  931bc7:	mov    rsi,rax
  931bca:	movq   xmm0,rax
  931bcf:	lea    rdx,[rcx+rcx*2]
  931bd3:	shl    rdx,0x3
  931bd7:	movdqu xmm1,XMMWORD PTR [r13+rdx*1-0x360]
  931be1:	sub    rsi,QWORD PTR [r12+rdx*1]
  931be5:	movzx  r13d,BYTE PTR [r12+0x181]
  931bee:	movq   xmm2,rsi
  931bf3:	punpcklqdq xmm0,xmm2
  931bf7:	psubq  xmm0,xmm1
  931bfb:	cmp    r13b,bl
  931bfe:	jae    931c20 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc10>
  931c00:	lea    eax,[r15-0x1]
  931c04:	movdqa xmm1,xmm0
  931c08:	cdqe
  931c0a:	lea    rax,[rax+rax*2]
  931c0e:	movq   xmm2,QWORD PTR [r12+rax*8+0x8]
  931c15:	paddq  xmm1,xmm2
  931c19:	movq   QWORD PTR [r12+rax*8+0x8],xmm1
  931c20:	movzx  r12d,WORD PTR [r14+0xfe0]
  931c28:	cmp    r12w,0x7f
  931c2d:	je     931da0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xd90>
  931c33:	mov    QWORD PTR [rbp-0x628],rcx
  931c3a:	sub    ebx,r13d
  931c3d:	movaps XMMWORD PTR [rbp-0x620],xmm0
  931c44:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  931c50:	mov    rcx,QWORD PTR [rbp-0x628]
  931c57:	movdqa xmm0,XMMWORD PTR [rbp-0x620]
  931c5f:	mov    DWORD PTR [rbp-0x1e8],0x0
  931c69:	lea    rdx,[rcx+rcx*2]
  931c6d:	movups XMMWORD PTR [rbp-0x1f8],xmm0
  931c74:	lea    rdx,[rax+rdx*8-0x360]
  931c7c:	mov    rax,QWORD PTR [rdx]
  931c7f:	mov    QWORD PTR [rbp-0x200],rax
  931c86:	movzx  eax,BYTE PTR [rdx+0x10]
  931c8a:	lea    edx,[r12+0x1]
  931c8f:	movdqa xmm2,XMMWORD PTR [rbp-0x200]
  931c97:	mov    WORD PTR [r14+0xfe0],dx
  931c9f:	mov    ah,bl
  931ca1:	mov    WORD PTR [rbp-0x1e4],ax
  931ca8:	movzx  eax,r12w
  931cac:	shl    rax,0x5
  931cb0:	add    rax,r14
  931cb3:	movups XMMWORD PTR [rax],xmm2
  931cb6:	movdqu xmm2,XMMWORD PTR [rbp-0x1f2]
  931cbe:	movups XMMWORD PTR [rax+0xe],xmm2
  931cc2:	jmp    931b93 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb83>
  931cc7:	shl    r8,0x4
  931ccb:	add    rdi,r8
  931cce:	cmp    r9,rdi
  931cd1:	je     931104 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xf4>
  931cd7:	mov    rcx,QWORD PTR [rbp-0x698]
  931cde:	mov    QWORD PTR [rcx+0x8],rdi
  931ce2:	jmp    931104 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xf4>
  931ce7:	mov    eax,0x1
  931cec:	lock xadd DWORD PTR [rip+0x2e1e6c],eax        # c13b60 <cpl::Profiling::registerRegion(char const*)::counter>
  931cf4:	mov    edx,0x1
  931cf9:	mov    r13d,0x1
  931cff:	add    eax,0x2
  931d02:	cmp    eax,0xfe
  931d07:	ja     931d22 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xd12>
  931d09:	lea    rdx,[rip+0x2e0310]        # c12020 <cpl::Profiling::regions>
  931d10:	mov    ecx,eax
  931d12:	lea    rdi,[rip+0xb642d]        # 9e8146 <_IO_stdin_used+0x18146>
  931d19:	mov    r13d,eax
  931d1c:	mov    QWORD PTR [rdx+rcx*8],rdi
  931d20:	mov    edx,eax
  931d22:	xor    eax,eax
  931d24:	lock cmpxchg BYTE PTR [rip+0x2ed3dc],dl        # c1f108 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  931d2c:	cmovne r13d,eax
  931d30:	jmp    931073 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x63>
  931d35:	mov    rdi,QWORD PTR [rbp-0x698]
  931d3c:	sub    rsi,rcx
  931d3f:	call   91d820 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  931d44:	mov    rax,QWORD PTR [rbp-0x6a8]
  931d4b:	mov    rdx,QWORD PTR [rax+0x38]
  931d4f:	mov    rax,QWORD PTR [rax+0x48]
  931d53:	jmp    931104 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xf4>
  931d58:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  931d5d:	mov    BYTE PTR [rbp-0x1f0],r13b
  931d64:	mov    QWORD PTR [rbp-0x200],rax
  931d6b:	lea    rax,[r12+r12*2]
  931d6f:	mov    QWORD PTR [rbp-0x1f8],0x0
  931d7a:	shl    rax,0x3
  931d7e:	movdqa xmm2,XMMWORD PTR [rbp-0x200]
  931d86:	movups XMMWORD PTR [rbx+rax*1-0x360],xmm2
  931d8e:	mov    BYTE PTR [rax+r14*1+0x10],r13b
  931d93:	movzx  r12d,BYTE PTR [r14+0x180]
  931d9b:	jmp    9310a5 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x95>
  931da0:	movzx  eax,WORD PTR [r14+0xfe2]
  931da8:	cmp    ax,0xffff
  931dac:	je     931b93 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb83>
  931db2:	add    eax,0x1
  931db5:	mov    WORD PTR [r14+0xfe2],ax
  931dbd:	jmp    931b93 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb83>
  931dc2:	call   920350 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  931dc7:	mov    rax,QWORD PTR [rbp-0x38]
  931dcb:	sub    rax,QWORD PTR fs:0x28
  931dd4:	je     931de4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xdd4>
  931dd6:	call   e0f90 <__stack_chk_fail@plt>
  931ddb:	endbr64
  931ddf:	mov    rbx,rax
  931de2:	jmp    931dc2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xdb2>
  931de4:	mov    rdi,rbx
  931de7:	call   e1b80 <_Unwind_Resume@plt>

Disassembly of section .fini:
