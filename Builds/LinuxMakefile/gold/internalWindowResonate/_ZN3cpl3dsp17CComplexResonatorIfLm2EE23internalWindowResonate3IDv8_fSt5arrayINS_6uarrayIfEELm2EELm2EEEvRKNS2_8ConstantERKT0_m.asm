; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000937e40 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long):
  937e40:	endbr64
  937e44:	lea    r10,[rsp+0x8]
  937e49:	and    rsp,0xffffffffffffffe0
  937e4d:	push   QWORD PTR [r10-0x8]
  937e51:	push   rbp
  937e52:	mov    rbp,rsp
  937e55:	push   r15
  937e57:	push   r14
  937e59:	push   r13
  937e5b:	push   r12
  937e5d:	push   r10
  937e5f:	push   rbx
  937e60:	sub    rsp,0x740
  937e67:	mov    QWORD PTR [rbp-0x698],rdi
  937e6e:	mov    QWORD PTR [rbp-0x6a8],rsi
  937e75:	mov    QWORD PTR [rbp-0x6b0],rdx
  937e7c:	mov    QWORD PTR [rbp-0x650],rcx
  937e83:	mov    rax,QWORD PTR fs:0x28
  937e8c:	mov    QWORD PTR [rbp-0x38],rax
  937e90:	xor    eax,eax
  937e92:	movzx  r13d,BYTE PTR [rip+0x2e7264]        # c1f0fe <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  937e9a:	test   r13b,r13b
  937e9d:	je     938b17 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xcd7>
  937ea3:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  937eaf:	cmp    QWORD PTR [rax-0x1d8],0x0
  937eb7:	mov    rbx,rax
  937eba:	lea    r14,[rax-0x360]
  937ec1:	je     937eec <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xac>
  937ec3:	movzx  r12d,BYTE PTR [r14+0x180]
  937ecb:	cmp    r12b,0xf
  937ecf:	jbe    938b88 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xd48>
  937ed5:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  937ee1:	add    r12d,0x1
  937ee5:	mov    BYTE PTR [rax-0x1e0],r12b
  937eec:	mov    rdx,QWORD PTR [rbp-0x6a8]
  937ef3:	mov    rcx,QWORD PTR [rbp-0x698]
  937efa:	mov    rax,QWORD PTR [rdx+0x48]
  937efe:	mov    rdx,QWORD PTR [rdx+0x38]
  937f02:	mov    r9,QWORD PTR [rcx+0x8]
  937f06:	mov    rdi,QWORD PTR [rcx]
  937f09:	mov    r8,rdx
  937f0c:	imul   r8,rax
  937f10:	mov    rcx,r9
  937f13:	sub    rcx,rdi
  937f16:	sar    rcx,0x2
  937f1a:	lea    rsi,[r8*4+0x0]
  937f22:	cmp    rcx,rsi
  937f25:	jb     938b65 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xd25>
  937f2b:	cmp    rsi,rcx
  937f2e:	jb     938af7 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xcb7>
  937f34:	mov    rdi,QWORD PTR [rbp-0x6a8]
  937f3b:	lea    rcx,[rax+rax*1]
  937f3f:	imul   rdx,rcx
  937f43:	cmp    QWORD PTR [rdi+0x40],0x0
  937f48:	je     93897e <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb3e>
  937f4e:	lea    rdi,[rcx+rax*1]
  937f52:	lea    r10,[rax*4+0x0]
  937f5a:	mov    QWORD PTR [rbp-0x680],0x0
  937f65:	shl    rdi,0x2
  937f69:	lea    rsi,[r10+rdx*1]
  937f6d:	mov    QWORD PTR [rbp-0x678],r10
  937f74:	lea    r13,[rbp-0xd0]
  937f7b:	mov    QWORD PTR [rbp-0x670],rdi
  937f82:	mov    rdi,rax
  937f85:	add    rax,r10
  937f88:	lea    r9,[rbp-0x1d0]
  937f8f:	shl    rax,0x2
  937f93:	shl    rdi,0x4
  937f97:	mov    QWORD PTR [rbp-0x6a0],0x0
  937fa2:	lea    r8,[rbp-0x190]
  937fa9:	mov    QWORD PTR [rbp-0x660],rax
  937fb0:	lea    rax,[rdx*4+0x0]
  937fb8:	lea    r12,[rbp-0x90]
  937fbf:	mov    QWORD PTR [rbp-0x690],rax
  937fc6:	mov    rax,rsi
  937fc9:	lea    rdx,[rbp-0x110]
  937fd0:	sub    rax,rcx
  937fd3:	mov    QWORD PTR [rbp-0x668],rdi
  937fda:	lea    rcx,[rbp-0x150]
  937fe1:	shl    rax,0x2
  937fe5:	mov    QWORD PTR [rbp-0x758],r10
  937fec:	mov    r14,rcx
  937fef:	mov    QWORD PTR [rbp-0x6c0],rax
  937ff6:	lea    rax,[rsi*4+0x0]
  937ffe:	mov    QWORD PTR [rbp-0x6b8],rax
  938005:	lea    rax,[rbp-0x470]
  93800c:	mov    QWORD PTR [rbp-0x6d8],rax
  938013:	lea    rax,[rbp-0x450]
  93801a:	mov    QWORD PTR [rbp-0x6f0],rax
  938021:	lea    rax,[rbp-0x430]
  938028:	mov    QWORD PTR [rbp-0x708],rax
  93802f:	lea    rax,[rbp-0x410]
  938036:	mov    QWORD PTR [rbp-0x720],rax
  93803d:	lea    rax,[rbp-0x3f0]
  938044:	mov    QWORD PTR [rbp-0x738],rax
  93804b:	lea    rax,[rbp-0x3d0]
  938052:	mov    QWORD PTR [rbp-0x750],rax
  938059:	lea    rax,[rbp-0x390]
  938060:	mov    QWORD PTR [rbp-0x6e8],rax
  938067:	lea    rax,[rbp-0x370]
  93806e:	mov    QWORD PTR [rbp-0x700],rax
  938075:	lea    rax,[rbp-0x350]
  93807c:	mov    QWORD PTR [rbp-0x718],rax
  938083:	lea    rax,[rbp-0x330]
  93808a:	mov    QWORD PTR [rbp-0x730],rax
  938091:	lea    rax,[rbp-0x310]
  938098:	mov    QWORD PTR [rbp-0x748],rax
  93809f:	lea    rax,[rbp-0x2f0]
  9380a6:	mov    QWORD PTR [rbp-0x6d0],rax
  9380ad:	lea    rax,[rbp-0x2d0]
  9380b4:	mov    QWORD PTR [rbp-0x6e0],rax
  9380bb:	lea    rax,[rbp-0x2b0]
  9380c2:	mov    QWORD PTR [rbp-0x6f8],rax
  9380c9:	lea    rax,[rbp-0x290]
  9380d0:	mov    QWORD PTR [rbp-0x710],rax
  9380d7:	lea    rax,[rbp-0x270]
  9380de:	mov    QWORD PTR [rbp-0x728],rax
  9380e5:	lea    rax,[rbp-0x250]
  9380ec:	mov    QWORD PTR [rbp-0x740],rax
  9380f3:	lea    rax,[rbp-0x230]
  9380fa:	mov    QWORD PTR [rbp-0x6c8],rax
  938101:	lea    rax,[rbp-0x3b0]
  938108:	mov    QWORD PTR [rbp-0x640],rax
  93810f:	mov    rax,r13
  938112:	nop    WORD PTR [rax+rax*1+0x0]
  938118:	mov    QWORD PTR [rbp-0x630],rax
  93811f:	mov    rax,QWORD PTR [rbp-0x6a8]
  938126:	mov    r13,QWORD PTR [rbp-0x680]
  93812d:	mov    rdi,QWORD PTR [rbp-0x6d8]
  938134:	mov    QWORD PTR [rbp-0x638],r8
  93813b:	mov    rbx,QWORD PTR [rax]
  93813e:	mov    QWORD PTR [rbp-0x628],r9
  938145:	mov    QWORD PTR [rbp-0x620],rdx
  93814c:	lea    rsi,[rbx+r13*1]
  938150:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  938155:	mov    r15,QWORD PTR [rbp-0x678]
  93815c:	mov    rdi,QWORD PTR [rbp-0x6f0]
  938163:	lea    rsi,[rbx+r15*1]
  938167:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  93816c:	mov    rax,QWORD PTR [rbp-0x758]
  938173:	mov    rdi,QWORD PTR [rbp-0x708]
  93817a:	add    rax,r15
  93817d:	lea    rsi,[rbx+rax*1]
  938181:	mov    QWORD PTR [rbp-0x688],rax
  938188:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  93818d:	mov    rdx,QWORD PTR [rbp-0x670]
  938194:	mov    rdi,QWORD PTR [rbp-0x720]
  93819b:	lea    rsi,[rbx+rdx*1]
  93819f:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9381a4:	mov    rcx,QWORD PTR [rbp-0x668]
  9381ab:	mov    rdi,QWORD PTR [rbp-0x738]
  9381b2:	lea    rsi,[rbx+rcx*1]
  9381b6:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9381bb:	mov    rdi,QWORD PTR [rbp-0x660]
  9381c2:	lea    rsi,[rbx+rdi*1]
  9381c6:	mov    rdi,QWORD PTR [rbp-0x750]
  9381cd:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9381d2:	mov    rbx,QWORD PTR [rbp-0x6b0]
  9381d9:	mov    rdi,QWORD PTR [rbp-0x6e8]
  9381e0:	mov    rsi,QWORD PTR [rbx]
  9381e3:	mov    QWORD PTR [rbp-0x1e0],rsi
  9381ea:	mov    rsi,QWORD PTR [rbp-0x698]
  9381f1:	mov    rbx,QWORD PTR [rsi]
  9381f4:	lea    rsi,[rbx+r13*1]
  9381f8:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9381fd:	mov    rdi,QWORD PTR [rbp-0x700]
  938204:	lea    rsi,[rbx+r15*1]
  938208:	movdqa xmm4,XMMWORD PTR [rbp-0x390]
  938210:	movdqa xmm3,XMMWORD PTR [rbp-0x380]
  938218:	movaps XMMWORD PTR [rbp-0x1d0],xmm4
  93821f:	movaps XMMWORD PTR [rbp-0x1c0],xmm3
  938226:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  93822b:	mov    rax,QWORD PTR [rbp-0x688]
  938232:	movdqa xmm7,XMMWORD PTR [rbp-0x370]
  93823a:	movdqa xmm2,XMMWORD PTR [rbp-0x360]
  938242:	mov    rdi,QWORD PTR [rbp-0x718]
  938249:	lea    rsi,[rbx+rax*1]
  93824d:	movaps XMMWORD PTR [rbp-0x190],xmm7
  938254:	movaps XMMWORD PTR [rbp-0x180],xmm2
  93825b:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  938260:	mov    rdx,QWORD PTR [rbp-0x670]
  938267:	movdqa xmm4,XMMWORD PTR [rbp-0x350]
  93826f:	movdqa xmm3,XMMWORD PTR [rbp-0x340]
  938277:	mov    rdi,QWORD PTR [rbp-0x730]
  93827e:	lea    rsi,[rbx+rdx*1]
  938282:	movaps XMMWORD PTR [rbp-0x150],xmm4
  938289:	movaps XMMWORD PTR [rbp-0x140],xmm3
  938290:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  938295:	mov    rcx,QWORD PTR [rbp-0x668]
  93829c:	movdqa xmm7,XMMWORD PTR [rbp-0x330]
  9382a4:	movdqa xmm2,XMMWORD PTR [rbp-0x320]
  9382ac:	mov    rdi,QWORD PTR [rbp-0x748]
  9382b3:	lea    rsi,[rbx+rcx*1]
  9382b7:	movaps XMMWORD PTR [rbp-0x110],xmm7
  9382be:	movaps XMMWORD PTR [rbp-0x100],xmm2
  9382c5:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9382ca:	mov    rdi,QWORD PTR [rbp-0x660]
  9382d1:	movdqa xmm4,XMMWORD PTR [rbp-0x310]
  9382d9:	movdqa xmm3,XMMWORD PTR [rbp-0x300]
  9382e1:	lea    rsi,[rbx+rdi*1]
  9382e5:	mov    rdi,QWORD PTR [rbp-0x6d0]
  9382ec:	movaps XMMWORD PTR [rbp-0xd0],xmm4
  9382f3:	movaps XMMWORD PTR [rbp-0xc0],xmm3
  9382fa:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9382ff:	mov    rax,QWORD PTR [rbp-0x6b0]
  938306:	movdqa xmm7,XMMWORD PTR [rbp-0x2f0]
  93830e:	movdqa xmm2,XMMWORD PTR [rbp-0x2e0]
  938316:	mov    rdi,QWORD PTR [rbp-0x6e0]
  93831d:	mov    rsi,QWORD PTR [rax+0x10]
  938321:	mov    rax,QWORD PTR [rbp-0x690]
  938328:	movaps XMMWORD PTR [rbp-0x90],xmm7
  93832f:	movaps XMMWORD PTR [rbp-0x80],xmm2
  938333:	mov    QWORD PTR [rbp-0x1d8],rsi
  93833a:	lea    rsi,[rax+r13*1]
  93833e:	add    rsi,rbx
  938341:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  938346:	mov    rax,QWORD PTR [rbp-0x690]
  93834d:	movdqa xmm4,XMMWORD PTR [rbp-0x2d0]
  938355:	movdqa xmm3,XMMWORD PTR [rbp-0x2c0]
  93835d:	mov    rdi,QWORD PTR [rbp-0x6f8]
  938364:	lea    rsi,[rax+r15*1]
  938368:	movaps XMMWORD PTR [rbp-0x1b0],xmm4
  93836f:	add    rsi,rbx
  938372:	movaps XMMWORD PTR [rbp-0x1a0],xmm3
  938379:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  93837e:	mov    rax,QWORD PTR [rbp-0x6c0]
  938385:	movdqa xmm7,XMMWORD PTR [rbp-0x2b0]
  93838d:	movdqa xmm2,XMMWORD PTR [rbp-0x2a0]
  938395:	mov    rdi,QWORD PTR [rbp-0x710]
  93839c:	lea    rsi,[rax+r13*1]
  9383a0:	movaps XMMWORD PTR [rbp-0x170],xmm7
  9383a7:	add    rsi,rbx
  9383aa:	movaps XMMWORD PTR [rbp-0x160],xmm2
  9383b1:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9383b6:	mov    rax,QWORD PTR [rbp-0x6c0]
  9383bd:	movdqa xmm3,XMMWORD PTR [rbp-0x280]
  9383c5:	movdqa xmm4,XMMWORD PTR [rbp-0x290]
  9383cd:	mov    rdi,QWORD PTR [rbp-0x728]
  9383d4:	lea    rsi,[rax+r15*1]
  9383d8:	movaps XMMWORD PTR [rbp-0x120],xmm3
  9383df:	add    rsi,rbx
  9383e2:	movaps XMMWORD PTR [rbp-0x130],xmm4
  9383e9:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  9383ee:	mov    rax,QWORD PTR [rbp-0x6b8]
  9383f5:	movdqa xmm7,XMMWORD PTR [rbp-0x270]
  9383fd:	movdqa xmm2,XMMWORD PTR [rbp-0x260]
  938405:	mov    rdi,QWORD PTR [rbp-0x740]
  93840c:	lea    rsi,[rax+r13*1]
  938410:	movaps XMMWORD PTR [rbp-0xf0],xmm7
  938417:	add    rsi,rbx
  93841a:	movaps XMMWORD PTR [rbp-0xe0],xmm2
  938421:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  938426:	mov    rax,QWORD PTR [rbp-0x6b8]
  93842d:	movdqa xmm5,XMMWORD PTR [rbp-0x250]
  938435:	movdqa xmm4,XMMWORD PTR [rbp-0x240]
  93843d:	mov    rdi,QWORD PTR [rbp-0x6c8]
  938444:	lea    rsi,[rax+r15*1]
  938448:	movaps XMMWORD PTR [rbp-0xb0],xmm5
  93844f:	add    rsi,rbx
  938452:	movaps XMMWORD PTR [rbp-0xa0],xmm4
  938459:	call   85fae0 <float __vector(8) cpl::simd::load<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  93845e:	movdqa xmm3,XMMWORD PTR [rbp-0x230]
  938466:	xor    esi,esi
  938468:	movdqa xmm7,XMMWORD PTR [rbp-0x220]
  938470:	cmp    QWORD PTR [rbp-0x650],0x0
  938478:	mov    rdx,QWORD PTR [rbp-0x620]
  93847f:	lea    rdi,[rbp-0x1e0]
  938486:	mov    r9,QWORD PTR [rbp-0x628]
  93848d:	mov    rax,QWORD PTR [rbp-0x630]
  938494:	movaps XMMWORD PTR [rbp-0x70],xmm3
  938498:	mov    r8,QWORD PTR [rbp-0x638]
  93849f:	movaps XMMWORD PTR [rbp-0x60],xmm7
  9384a3:	je     938796 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x956>
  9384a9:	mov    QWORD PTR [rbp-0x658],rdi
  9384b0:	mov    r15,r14
  9384b3:	mov    r14,r12
  9384b6:	mov    r12,rdx
  9384b9:	mov    rdx,rsi
  9384bc:	nop    DWORD PTR [rax+0x0]
  9384c0:	mov    rcx,QWORD PTR [rbp-0x658]
  9384c7:	mov    rdi,r15
  9384ca:	mov    QWORD PTR [rbp-0x648],rdx
  9384d1:	mov    rdx,r12
  9384d4:	mov    r15,r14
  9384d7:	xor    ebx,ebx
  9384d9:	mov    r14,rdi
  9384dc:	mov    r12,rcx
  9384df:	mov    rcx,rdx
  9384e2:	mov    r13,QWORD PTR [r12]
  9384e6:	mov    rdi,QWORD PTR [rbp-0x640]
  9384ed:	mov    QWORD PTR [rbp-0x638],r8
  9384f4:	add    r12,0x8
  9384f8:	mov    QWORD PTR [rbp-0x630],rax
  9384ff:	mov    rsi,r13
  938502:	mov    QWORD PTR [rbp-0x628],r9
  938509:	add    r13,0x4
  93850d:	mov    QWORD PTR [rbp-0x620],rcx
  938514:	call   9243d0 <float __vector(8) cpl::simd::broadcast<float __vector(8)>(cpl::simd::scalar_of<float __vector(8), 4ul>::type const*)>
  938519:	mov    r9,QWORD PTR [rbp-0x628]
  938520:	mov    r8,QWORD PTR [rbp-0x638]
  938527:	movaps xmm0,XMMWORD PTR [rbp-0x450]
  93852e:	movaps xmm2,XMMWORD PTR [rbp-0x470]
  938535:	lea    rsi,[r9+rbx*1]
  938539:	movaps xmm3,XMMWORD PTR [rbp-0x460]
  938540:	mov    rcx,QWORD PTR [rbp-0x620]
  938547:	movaps xmm1,XMMWORD PTR [rsi+0x10]
  93854b:	movaps xmm6,XMMWORD PTR [rsi]
  93854e:	lea    rsi,[r8+rbx*1]
  938552:	movaps xmm5,XMMWORD PTR [rsi]
  938555:	movaps xmm4,XMMWORD PTR [rsi+0x10]
  938559:	lea    rsi,[r14+rbx*1]
  93855d:	mulps  xmm0,xmm6
  938560:	mulps  xmm2,xmm5
  938563:	mulps  xmm3,xmm4
  938566:	mulps  xmm6,XMMWORD PTR [rbp-0x470]
  93856d:	mulps  xmm5,XMMWORD PTR [rbp-0x450]
  938574:	mulps  xmm4,XMMWORD PTR [rbp-0x440]
  93857b:	addps  xmm0,xmm2
  93857e:	movaps xmm2,XMMWORD PTR [rbp-0x440]
  938585:	mulps  xmm2,xmm1
  938588:	mulps  xmm1,XMMWORD PTR [rbp-0x460]
  93858f:	addps  xmm2,xmm3
  938592:	subps  xmm1,xmm4
  938595:	movaps XMMWORD PTR [rbp-0x600],xmm2
  93859c:	movaps XMMWORD PTR [r8+rbx*1],xmm0
  9385a1:	movaps xmm0,xmm6
  9385a4:	movdqa xmm7,XMMWORD PTR [rbp-0x600]
  9385ac:	subps  xmm0,xmm5
  9385af:	movaps XMMWORD PTR [r8+rbx*1+0x10],xmm7
  9385b5:	movaps xmm3,XMMWORD PTR [rbp-0x3b0]
  9385bc:	movaps xmm2,XMMWORD PTR [rbp-0x3a0]
  9385c3:	addps  xmm0,xmm3
  9385c6:	addps  xmm1,xmm2
  9385c9:	movaps XMMWORD PTR [rbp-0x610],xmm0
  9385d0:	movaps XMMWORD PTR [rbp-0x600],xmm1
  9385d7:	movaps xmm1,XMMWORD PTR [rbp-0x410]
  9385de:	movaps XMMWORD PTR [r9+rbx*1],xmm0
  9385e3:	movdqa xmm4,XMMWORD PTR [rbp-0x600]
  9385eb:	movaps xmm7,XMMWORD PTR [rsi]
  9385ee:	movaps xmm0,XMMWORD PTR [rbp-0x430]
  9385f5:	movaps XMMWORD PTR [r9+rbx*1+0x10],xmm4
  9385fb:	movaps xmm5,XMMWORD PTR [rsi+0x10]
  9385ff:	lea    rsi,[rcx+rbx*1]
  938603:	mulps  xmm1,xmm7
  938606:	movaps xmm6,XMMWORD PTR [rsi]
  938609:	movaps xmm4,XMMWORD PTR [rsi+0x10]
  93860d:	mulps  xmm0,xmm6
  938610:	addps  xmm0,xmm1
  938613:	movaps xmm1,XMMWORD PTR [rbp-0x420]
  93861a:	movaps xmm8,XMMWORD PTR [rbp-0x400]
  938622:	mulps  xmm7,XMMWORD PTR [rbp-0x430]
  938629:	mov    rax,QWORD PTR [rbp-0x630]
  938630:	mulps  xmm8,xmm5
  938634:	mulps  xmm1,xmm4
  938637:	lea    rsi,[rax+rbx*1]
  93863b:	mulps  xmm5,XMMWORD PTR [rbp-0x420]
  938642:	mulps  xmm6,XMMWORD PTR [rbp-0x410]
  938649:	mulps  xmm4,XMMWORD PTR [rbp-0x400]
  938650:	addps  xmm1,xmm8
  938654:	movaps xmm8,XMMWORD PTR [rbp-0x3c0]
  93865c:	movaps XMMWORD PTR [rbp-0x600],xmm1
  938663:	movaps xmm1,xmm5
  938666:	movaps XMMWORD PTR [rcx+rbx*1],xmm0
  93866a:	movdqa xmm0,XMMWORD PTR [rbp-0x600]
  938672:	subps  xmm1,xmm4
  938675:	movaps XMMWORD PTR [rcx+rbx*1+0x10],xmm0
  93867a:	movaps xmm0,xmm7
  93867d:	subps  xmm0,xmm6
  938680:	addps  xmm1,xmm2
  938683:	addps  xmm0,xmm3
  938686:	movaps XMMWORD PTR [rbp-0x600],xmm1
  93868d:	movaps xmm1,XMMWORD PTR [rbp-0x3d0]
  938694:	movaps XMMWORD PTR [rbp-0x610],xmm0
  93869b:	movaps XMMWORD PTR [r14+rbx*1],xmm0
  9386a0:	movdqa xmm7,XMMWORD PTR [rbp-0x600]
  9386a8:	movaps xmm0,XMMWORD PTR [rbp-0x3f0]
  9386af:	movaps XMMWORD PTR [r14+rbx*1+0x10],xmm7
  9386b5:	movaps xmm7,XMMWORD PTR [rsi]
  9386b8:	movaps xmm5,XMMWORD PTR [rsi+0x10]
  9386bc:	lea    rsi,[r15+rbx*1]
  9386c0:	movaps xmm6,XMMWORD PTR [rsi]
  9386c3:	movaps xmm4,XMMWORD PTR [rsi+0x10]
  9386c7:	mulps  xmm1,xmm7
  9386ca:	mulps  xmm0,xmm6
  9386cd:	mulps  xmm8,xmm5
  9386d1:	mulps  xmm7,XMMWORD PTR [rbp-0x3f0]
  9386d8:	mulps  xmm6,XMMWORD PTR [rbp-0x3d0]
  9386df:	mulps  xmm5,XMMWORD PTR [rbp-0x3e0]
  9386e6:	addps  xmm0,xmm1
  9386e9:	movaps xmm1,XMMWORD PTR [rbp-0x3e0]
  9386f0:	mulps  xmm1,xmm4
  9386f3:	addps  xmm1,xmm8
  9386f7:	movaps XMMWORD PTR [rbp-0x600],xmm1
  9386fe:	movaps xmm1,xmm5
  938701:	movaps XMMWORD PTR [r15+rbx*1],xmm0
  938706:	movdqa xmm0,XMMWORD PTR [rbp-0x600]
  93870e:	movaps XMMWORD PTR [r15+rbx*1+0x10],xmm0
  938714:	movaps xmm0,xmm7
  938717:	subps  xmm0,xmm6
  93871a:	addps  xmm0,xmm3
  93871d:	movaps XMMWORD PTR [rbp-0x610],xmm0
  938724:	movaps xmm3,XMMWORD PTR [rbp-0x3c0]
  93872b:	mulps  xmm3,xmm4
  93872e:	subps  xmm1,xmm3
  938731:	addps  xmm1,xmm2
  938734:	movaps XMMWORD PTR [rbp-0x600],xmm1
  93873b:	movaps XMMWORD PTR [rax+rbx*1],xmm0
  93873f:	movdqa xmm3,XMMWORD PTR [rbp-0x600]
  938747:	mov    QWORD PTR [r12-0x8],r13
  93874c:	movaps XMMWORD PTR [rax+rbx*1+0x10],xmm3
  938751:	add    rbx,0x20
  938755:	cmp    rbx,0x40
  938759:	jne    9384e2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x6a2>
  93875f:	mov    rdx,QWORD PTR [rbp-0x648]
  938766:	mov    rdi,r14
  938769:	mov    r12,rcx
  93876c:	mov    r14,r15
  93876f:	mov    r15,rdi
  938772:	add    rdx,0x1
  938776:	cmp    QWORD PTR [rbp-0x650],rdx
  93877d:	jne    9384c0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x680>
  938783:	mov    rdx,rcx
  938786:	mov    rcx,QWORD PTR [rbp-0x698]
  93878d:	mov    r12,r14
  938790:	mov    r14,rdi
  938793:	mov    rbx,QWORD PTR [rcx]
  938796:	xor    r15d,r15d
  938799:	mov    rdi,rbx
  93879c:	mov    rbx,r14
  93879f:	mov    r14,r12
  9387a2:	mov    r12,r15
  9387a5:	mov    r15,QWORD PTR [rbp-0x698]
  9387ac:	xor    r13d,r13d
  9387af:	movdqa xmm4,XMMWORD PTR [r9+r12*1]
  9387b5:	sub    rsp,0x20
  9387b9:	mov    QWORD PTR [rbp-0x648],rax
  9387c0:	mov    rax,QWORD PTR [rbp-0x680]
  9387c7:	mov    QWORD PTR [rbp-0x620],rdx
  9387ce:	mov    QWORD PTR [rbp-0x628],r8
  9387d5:	movaps XMMWORD PTR [rsp],xmm4
  9387d9:	movdqa xmm4,XMMWORD PTR [r9+r12*1+0x10]
  9387e0:	lea    rsi,[rax+r13*1]
  9387e4:	mov    QWORD PTR [rbp-0x638],r9
  9387eb:	add    rdi,rsi
  9387ee:	movaps XMMWORD PTR [rsp+0x10],xmm4
  9387f3:	call   85fb10 <cpl::simd::store(float*, float __vector(8))>
  9387f8:	mov    r8,QWORD PTR [rbp-0x628]
  9387ff:	mov    rax,QWORD PTR [rbp-0x678]
  938806:	movdqa xmm4,XMMWORD PTR [r8+r12*1]
  93880c:	lea    rdi,[rax+r13*1]
  938810:	add    rdi,QWORD PTR [r15]
  938813:	movaps XMMWORD PTR [rsp],xmm4
  938817:	movdqa xmm4,XMMWORD PTR [r8+r12*1+0x10]
  93881e:	mov    QWORD PTR [rbp-0x630],r8
  938825:	movaps XMMWORD PTR [rsp+0x10],xmm4
  93882a:	call   85fb10 <cpl::simd::store(float*, float __vector(8))>
  93882f:	movdqa xmm4,XMMWORD PTR [rbx+r12*1]
  938835:	mov    rax,QWORD PTR [rbp-0x688]
  93883c:	lea    rdi,[rax+r13*1]
  938840:	add    rdi,QWORD PTR [r15]
  938843:	movaps XMMWORD PTR [rsp],xmm4
  938847:	movdqa xmm3,XMMWORD PTR [rbx+r12*1+0x10]
  93884e:	movaps XMMWORD PTR [rsp+0x10],xmm3
  938853:	call   85fb10 <cpl::simd::store(float*, float __vector(8))>
  938858:	mov    rdx,QWORD PTR [rbp-0x620]
  93885f:	mov    rax,QWORD PTR [rbp-0x670]
  938866:	movdqa xmm7,XMMWORD PTR [rdx+r12*1]
  93886c:	lea    rdi,[rax+r13*1]
  938870:	add    rdi,QWORD PTR [r15]
  938873:	movaps XMMWORD PTR [rsp],xmm7
  938877:	movdqa xmm2,XMMWORD PTR [rdx+r12*1+0x10]
  93887e:	mov    QWORD PTR [rbp-0x628],rdx
  938885:	movaps XMMWORD PTR [rsp+0x10],xmm2
  93888a:	call   85fb10 <cpl::simd::store(float*, float __vector(8))>
  93888f:	mov    rax,QWORD PTR [rbp-0x668]
  938896:	lea    rdi,[rax+r13*1]
  93889a:	mov    rax,QWORD PTR [rbp-0x648]
  9388a1:	add    rdi,QWORD PTR [r15]
  9388a4:	movdqa xmm4,XMMWORD PTR [rax+r12*1]
  9388aa:	movaps XMMWORD PTR [rsp],xmm4
  9388ae:	movdqa xmm3,XMMWORD PTR [rax+r12*1+0x10]
  9388b5:	mov    QWORD PTR [rbp-0x620],rax
  9388bc:	movaps XMMWORD PTR [rsp+0x10],xmm3
  9388c1:	call   85fb10 <cpl::simd::store(float*, float __vector(8))>
  9388c6:	mov    rax,QWORD PTR [rbp-0x660]
  9388cd:	movdqa xmm7,XMMWORD PTR [r14+r12*1]
  9388d3:	lea    rdi,[rax+r13*1]
  9388d7:	add    rdi,QWORD PTR [r15]
  9388da:	movaps XMMWORD PTR [rsp],xmm7
  9388de:	movdqa xmm2,XMMWORD PTR [r14+r12*1+0x10]
  9388e5:	add    r12,0x20
  9388e9:	movaps XMMWORD PTR [rsp+0x10],xmm2
  9388ee:	call   85fb10 <cpl::simd::store(float*, float __vector(8))>
  9388f3:	add    rsp,0x20
  9388f7:	cmp    r12,0x40
  9388fb:	mov    rax,QWORD PTR [rbp-0x620]
  938902:	mov    rdx,QWORD PTR [rbp-0x628]
  938909:	mov    r8,QWORD PTR [rbp-0x630]
  938910:	mov    r9,QWORD PTR [rbp-0x638]
  938917:	je     938930 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xaf0>
  938919:	mov    rcx,QWORD PTR [rbp-0x690]
  938920:	mov    rdi,QWORD PTR [r15]
  938923:	add    r13,rcx
  938926:	jmp    9387af <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x96f>
  93892b:	nop    DWORD PTR [rax+rax*1+0x0]
  938930:	mov    rdi,QWORD PTR [rbp-0x6a8]
  938937:	add    QWORD PTR [rbp-0x6a0],0x8
  93893f:	mov    r12,r14
  938942:	mov    r14,rbx
  938945:	add    QWORD PTR [rbp-0x680],0x20
  93894d:	mov    rcx,QWORD PTR [rbp-0x6a0]
  938954:	add    QWORD PTR [rbp-0x678],0x20
  93895c:	add    QWORD PTR [rbp-0x670],0x20
  938964:	add    QWORD PTR [rbp-0x668],0x20
  93896c:	add    QWORD PTR [rbp-0x660],0x20
  938974:	cmp    rcx,QWORD PTR [rdi+0x40]
  938978:	jb     938118 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x2d8>
  93897e:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  93898a:	lea    r12,[rax-0x360]
  938991:	mov    r13,rax
  938994:	mov    r14,QWORD PTR [r12+0x188]
  93899c:	test   r14,r14
  93899f:	je     9389c3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb83>
  9389a1:	movzx  eax,BYTE PTR [r12+0x180]
  9389aa:	lea    ebx,[rax-0x1]
  9389ad:	mov    BYTE PTR [r12+0x180],bl
  9389b5:	cmp    bl,0xf
  9389b8:	jbe    9389eb <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xbab>
  9389ba:	add    WORD PTR [r14+0xfe2],0x1
  9389c3:	mov    rax,QWORD PTR [rbp-0x38]
  9389c7:	sub    rax,QWORD PTR fs:0x28
  9389d0:	jne    938c06 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xdc6>
  9389d6:	lea    rsp,[rbp-0x30]
  9389da:	pop    rbx
  9389db:	pop    r10
  9389dd:	pop    r12
  9389df:	pop    r13
  9389e1:	pop    r14
  9389e3:	pop    r15
  9389e5:	pop    rbp
  9389e6:	lea    rsp,[r10-0x8]
  9389ea:	ret
  9389eb:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  9389f0:	movzx  r15d,bl
  9389f4:	movsxd rcx,r15d
  9389f7:	mov    rsi,rax
  9389fa:	movq   xmm0,rax
  9389ff:	lea    rdx,[rcx+rcx*2]
  938a03:	shl    rdx,0x3
  938a07:	movdqu xmm1,XMMWORD PTR [r13+rdx*1-0x360]
  938a11:	sub    rsi,QWORD PTR [r12+rdx*1]
  938a15:	movzx  r13d,BYTE PTR [r12+0x181]
  938a1e:	movq   xmm2,rsi
  938a23:	punpcklqdq xmm0,xmm2
  938a27:	psubq  xmm0,xmm1
  938a2b:	cmp    r13b,bl
  938a2e:	jae    938a50 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc10>
  938a30:	lea    eax,[r15-0x1]
  938a34:	movdqa xmm1,xmm0
  938a38:	cdqe
  938a3a:	lea    rax,[rax+rax*2]
  938a3e:	movq   xmm2,QWORD PTR [r12+rax*8+0x8]
  938a45:	paddq  xmm1,xmm2
  938a49:	movq   QWORD PTR [r12+rax*8+0x8],xmm1
  938a50:	movzx  r12d,WORD PTR [r14+0xfe0]
  938a58:	cmp    r12w,0x7f
  938a5d:	je     938bd0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xd90>
  938a63:	mov    QWORD PTR [rbp-0x628],rcx
  938a6a:	sub    ebx,r13d
  938a6d:	movaps XMMWORD PTR [rbp-0x620],xmm0
  938a74:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  938a80:	mov    rcx,QWORD PTR [rbp-0x628]
  938a87:	movdqa xmm0,XMMWORD PTR [rbp-0x620]
  938a8f:	mov    DWORD PTR [rbp-0x1e8],0x0
  938a99:	lea    rdx,[rcx+rcx*2]
  938a9d:	movups XMMWORD PTR [rbp-0x1f8],xmm0
  938aa4:	lea    rdx,[rax+rdx*8-0x360]
  938aac:	mov    rax,QWORD PTR [rdx]
  938aaf:	mov    QWORD PTR [rbp-0x200],rax
  938ab6:	movzx  eax,BYTE PTR [rdx+0x10]
  938aba:	lea    edx,[r12+0x1]
  938abf:	movdqa xmm2,XMMWORD PTR [rbp-0x200]
  938ac7:	mov    WORD PTR [r14+0xfe0],dx
  938acf:	mov    ah,bl
  938ad1:	mov    WORD PTR [rbp-0x1e4],ax
  938ad8:	movzx  eax,r12w
  938adc:	shl    rax,0x5
  938ae0:	add    rax,r14
  938ae3:	movups XMMWORD PTR [rax],xmm2
  938ae6:	movdqu xmm2,XMMWORD PTR [rbp-0x1f2]
  938aee:	movups XMMWORD PTR [rax+0xe],xmm2
  938af2:	jmp    9389c3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb83>
  938af7:	shl    r8,0x4
  938afb:	add    rdi,r8
  938afe:	cmp    r9,rdi
  938b01:	je     937f34 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xf4>
  938b07:	mov    rcx,QWORD PTR [rbp-0x698]
  938b0e:	mov    QWORD PTR [rcx+0x8],rdi
  938b12:	jmp    937f34 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xf4>
  938b17:	mov    eax,0x1
  938b1c:	lock xadd DWORD PTR [rip+0x2db03c],eax        # c13b60 <cpl::Profiling::registerRegion(char const*)::counter>
  938b24:	mov    edx,0x1
  938b29:	mov    r13d,0x1
  938b2f:	add    eax,0x2
  938b32:	cmp    eax,0xfe
  938b37:	ja     938b52 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xd12>
  938b39:	lea    rdx,[rip+0x2d94e0]        # c12020 <cpl::Profiling::regions>
  938b40:	mov    ecx,eax
  938b42:	lea    rdi,[rip+0xaf5fd]        # 9e8146 <_IO_stdin_used+0x18146>
  938b49:	mov    r13d,eax
  938b4c:	mov    QWORD PTR [rdx+rcx*8],rdi
  938b50:	mov    edx,eax
  938b52:	xor    eax,eax
  938b54:	lock cmpxchg BYTE PTR [rip+0x2e65a2],dl        # c1f0fe <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  938b5c:	cmovne r13d,eax
  938b60:	jmp    937ea3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x63>
  938b65:	mov    rdi,QWORD PTR [rbp-0x698]
  938b6c:	sub    rsi,rcx
  938b6f:	call   91d820 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  938b74:	mov    rax,QWORD PTR [rbp-0x6a8]
  938b7b:	mov    rdx,QWORD PTR [rax+0x38]
  938b7f:	mov    rax,QWORD PTR [rax+0x48]
  938b83:	jmp    937f34 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xf4>
  938b88:	call   e25b0 <std::chrono::_V2::steady_clock::now()@plt>
  938b8d:	mov    BYTE PTR [rbp-0x1f0],r13b
  938b94:	mov    QWORD PTR [rbp-0x200],rax
  938b9b:	lea    rax,[r12+r12*2]
  938b9f:	mov    QWORD PTR [rbp-0x1f8],0x0
  938baa:	shl    rax,0x3
  938bae:	movdqa xmm2,XMMWORD PTR [rbp-0x200]
  938bb6:	movups XMMWORD PTR [rbx+rax*1-0x360],xmm2
  938bbe:	mov    BYTE PTR [rax+r14*1+0x10],r13b
  938bc3:	movzx  r12d,BYTE PTR [r14+0x180]
  938bcb:	jmp    937ed5 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x95>
  938bd0:	movzx  eax,WORD PTR [r14+0xfe2]
  938bd8:	cmp    ax,0xffff
  938bdc:	je     9389c3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb83>
  938be2:	add    eax,0x1
  938be5:	mov    WORD PTR [r14+0xfe2],ax
  938bed:	jmp    9389c3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xb83>
  938bf2:	call   920350 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  938bf7:	mov    rax,QWORD PTR [rbp-0x38]
  938bfb:	sub    rax,QWORD PTR fs:0x28
  938c04:	je     938c14 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xdd4>
  938c06:	call   e0f90 <__stack_chk_fail@plt>
  938c0b:	endbr64
  938c0f:	mov    rbx,rax
  938c12:	jmp    938bf2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xdb2>
  938c14:	mov    rdi,rbx
  938c17:	call   e1b80 <_Unwind_Resume@plt>

Disassembly of section .fini:
