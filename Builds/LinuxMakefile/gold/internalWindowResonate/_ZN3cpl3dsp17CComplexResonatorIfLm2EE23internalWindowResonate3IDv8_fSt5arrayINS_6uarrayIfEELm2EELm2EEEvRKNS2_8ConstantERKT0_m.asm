; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000953620 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long):
  953620:	endbr64
  953624:	push   rbp
  953625:	mov    rbp,rsp
  953628:	push   r15
  95362a:	push   r14
  95362c:	mov    r14,rdi
  95362f:	push   r13
  953631:	push   r12
  953633:	push   rbx
  953634:	mov    rbx,rcx
  953637:	and    rsp,0xffffffffffffffe0
  95363b:	sub    rsp,0x160
  953642:	mov    QWORD PTR [rsp+0x70],rsi
  953647:	mov    QWORD PTR [rsp+0x30],rdx
  95364c:	movzx  r13d,BYTE PTR [rip+0x2deaaa]        # c320fe <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  953654:	test   r13b,r13b
  953657:	je     953c28 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x608>
  95365d:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  953669:	cmp    QWORD PTR [rax-0x1d8],0x0
  953671:	lea    r15,[rax-0x360]
  953678:	je     9536a3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x83>
  95367a:	movzx  r12d,BYTE PTR [r15+0x180]
  953682:	cmp    r12b,0xf
  953686:	jbe    953c93 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x673>
  95368c:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  953698:	add    r12d,0x1
  95369c:	mov    BYTE PTR [rax-0x1e0],r12b
  9536a3:	mov    rax,QWORD PTR [rsp+0x70]
  9536a8:	mov    r9,QWORD PTR [r14+0x8]
  9536ac:	mov    rdi,QWORD PTR [r14]
  9536af:	mov    rdx,QWORD PTR [rax+0x48]
  9536b3:	mov    rax,QWORD PTR [rax+0x38]
  9536b7:	mov    rcx,r9
  9536ba:	sub    rcx,rdi
  9536bd:	sar    rcx,0x2
  9536c1:	mov    r8,rax
  9536c4:	imul   r8,rdx
  9536c8:	lea    rsi,[r8*4+0x0]
  9536d0:	cmp    rcx,rsi
  9536d3:	jb     953c76 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x656>
  9536d9:	cmp    rsi,rcx
  9536dc:	jb     953c0f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x5ef>
  9536e2:	mov    rdi,QWORD PTR [rsp+0x70]
  9536e7:	lea    rcx,[rdx+rdx*1]
  9536eb:	imul   rax,rcx
  9536ef:	cmp    QWORD PTR [rdi+0x40],0x0
  9536f4:	je     953aaa <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x48a>
  9536fa:	lea    r11,[rdx*4+0x0]
  953702:	xor    r10d,r10d
  953705:	mov    QWORD PTR [rsp+0x78],0x0
  95370e:	lea    rsi,[r11+rax*1]
  953712:	shl    rax,0x2
  953716:	mov    QWORD PTR [rsp+0x28],r11
  95371b:	mov    QWORD PTR [rsp+0x10],rax
  953720:	lea    rdi,[rdx*8+0x0]
  953728:	mov    rax,rsi
  95372b:	shl    rdx,0x4
  95372f:	sub    rax,rcx
  953732:	mov    QWORD PTR [rsp+0x20],rdi
  953737:	shl    rax,0x2
  95373b:	mov    QWORD PTR [rsp+0x18],rdx
  953740:	mov    QWORD PTR [rsp+0x60],rax
  953745:	lea    rax,[rsi*4+0x0]
  95374d:	mov    QWORD PTR [rsp+0x68],rax
  953752:	nop    WORD PTR [rax+rax*1+0x0]
  953758:	mov    rdi,QWORD PTR [rsp+0x20]
  95375d:	mov    rcx,QWORD PTR [rsp+0x28]
  953762:	mov    rax,QWORD PTR [rsp+0x70]
  953767:	mov    rdx,QWORD PTR [rsp+0x30]
  95376c:	add    rdi,r10
  95376f:	mov    r12,QWORD PTR [rsp+0x10]
  953774:	lea    rsi,[rcx+rdi*1]
  953778:	lea    r9,[rcx+rsi*1]
  95377c:	mov    rcx,QWORD PTR [rsp+0x18]
  953781:	mov    rax,QWORD PTR [rax]
  953784:	lea    r13,[r10+r12*1]
  953788:	add    r12,r11
  95378b:	lea    r8,[r11+rcx*1]
  95378f:	mov    rcx,QWORD PTR [r14]
  953792:	vmovaps ymm6,YMMWORD PTR [rax+r10*1]
  953798:	vmovaps ymm7,YMMWORD PTR [rax+r11*1]
  95379e:	vmovaps ymm9,YMMWORD PTR [rax+rdi*1]
  9537a3:	vmovaps ymm10,YMMWORD PTR [rax+rsi*1]
  9537a8:	lea    r15,[rcx+r10*1]
  9537ac:	vmovaps ymm11,YMMWORD PTR [rax+r9*1]
  9537b2:	vmovaps ymm8,YMMWORD PTR [rax+r8*1]
  9537b8:	mov    rax,QWORD PTR [rdx]
  9537bb:	vmovaps ymm2,YMMWORD PTR [rcx+r11*1]
  9537c1:	vmovaps ymm3,YMMWORD PTR [rcx+rsi*1]
  9537c6:	vmovaps YMMWORD PTR [rsp+0xa0],ymm2
  9537cf:	vmovaps ymm2,YMMWORD PTR [rcx+r9*1]
  9537d5:	mov    QWORD PTR [rsp+0x38],rax
  9537da:	mov    rax,QWORD PTR [rsp+0x60]
  9537df:	vmovaps YMMWORD PTR [rsp+0x100],ymm3
  9537e8:	vmovaps ymm3,YMMWORD PTR [rcx+r8*1]
  9537ee:	vmovaps YMMWORD PTR [rsp+0x120],ymm2
  9537f7:	vmovaps ymm2,YMMWORD PTR [rcx+r12*1]
  9537fd:	vmovaps YMMWORD PTR [rsp+0xc0],ymm3
  953806:	mov    rdx,QWORD PTR [rdx+0x10]
  95380a:	add    rax,r10
  95380d:	vmovaps ymm1,YMMWORD PTR [r15]
  953812:	vmovaps YMMWORD PTR [rsp+0x80],ymm2
  95381b:	mov    QWORD PTR [rsp+0x58],rax
  953820:	vmovaps ymm0,YMMWORD PTR [rcx+rdi*1]
  953825:	vmovaps ymm5,YMMWORD PTR [rcx+r13*1]
  95382b:	vmovaps ymm4,YMMWORD PTR [rcx+rax*1]
  953830:	mov    rax,QWORD PTR [rsp+0x60]
  953835:	add    rax,r11
  953838:	vmovaps ymm3,YMMWORD PTR [rcx+rax*1]
  95383d:	mov    QWORD PTR [rsp+0x50],rax
  953842:	mov    rax,QWORD PTR [rsp+0x68]
  953847:	vmovaps YMMWORD PTR [rsp+0xe0],ymm3
  953850:	add    rax,r10
  953853:	vmovaps ymm3,YMMWORD PTR [rcx+rax*1]
  953858:	mov    QWORD PTR [rsp+0x48],rax
  95385d:	mov    rax,QWORD PTR [rsp+0x68]
  953862:	add    rax,r11
  953865:	vmovaps ymm2,YMMWORD PTR [rcx+rax*1]
  95386a:	mov    QWORD PTR [rsp+0x40],rax
  95386f:	test   rbx,rbx
  953872:	je     9539d9 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3b9>
  953878:	mov    rax,QWORD PTR [rsp+0x38]
  95387d:	xor    r15d,r15d
  953880:	vmulps ymm15,ymm1,ymm7
  953884:	vmulps ymm1,ymm1,ymm6
  953888:	add    r15,0x1
  95388c:	add    rdx,0x4
  953890:	vbroadcastss ymm12,DWORD PTR [rax]
  953895:	add    rax,0x4
  953899:	vmovaps ymm14,YMMWORD PTR [rsp+0xa0]
  9538a2:	vaddps ymm1,ymm1,ymm12
  9538a7:	vmulps ymm13,ymm7,ymm14
  9538ac:	vmulps ymm14,ymm6,ymm14
  9538b1:	vaddps ymm15,ymm15,ymm14
  9538b6:	vmovaps ymm14,YMMWORD PTR [rsp+0x100]
  9538bf:	vsubps ymm1,ymm1,ymm13
  9538c4:	vmulps ymm13,ymm10,ymm14
  9538c9:	vmulps ymm14,ymm9,ymm14
  9538ce:	vmovaps YMMWORD PTR [rsp+0xa0],ymm15
  9538d7:	vmulps ymm15,ymm10,ymm0
  9538db:	vmulps ymm0,ymm9,ymm0
  9538df:	vaddps ymm15,ymm15,ymm14
  9538e4:	vmovaps ymm14,YMMWORD PTR [rsp+0xc0]
  9538ed:	vaddps ymm0,ymm0,ymm12
  9538f2:	vmovaps YMMWORD PTR [rsp+0x100],ymm15
  9538fb:	vmulps ymm15,ymm8,YMMWORD PTR [rsp+0x120]
  953904:	vsubps ymm0,ymm0,ymm13
  953909:	vmulps ymm13,ymm8,ymm14
  95390e:	vmulps ymm14,ymm11,ymm14
  953913:	vaddps ymm14,ymm15,ymm14
  953918:	vmovaps ymm15,YMMWORD PTR [rsp+0x80]
  953921:	vmovaps YMMWORD PTR [rsp+0xc0],ymm14
  95392a:	vmulps ymm14,ymm11,YMMWORD PTR [rsp+0x120]
  953933:	vaddps ymm14,ymm14,ymm12
  953938:	vmulps ymm12,ymm7,ymm15
  95393d:	vsubps ymm13,ymm14,ymm13
  953942:	vmulps ymm14,ymm6,ymm5
  953946:	vmulps ymm5,ymm7,ymm5
  95394a:	vsubps ymm12,ymm14,ymm12
  95394f:	vmulps ymm14,ymm6,ymm15
  953954:	vmovaps ymm15,YMMWORD PTR [rsp+0xe0]
  95395d:	vaddps ymm5,ymm5,ymm14
  953962:	vmovaps YMMWORD PTR [rsp+0x120],ymm13
  95396b:	vbroadcastss ymm13,DWORD PTR [rdx-0x4]
  953971:	vmulps ymm14,ymm10,ymm15
  953976:	vmovaps YMMWORD PTR [rsp+0x80],ymm5
  95397f:	vaddps ymm5,ymm13,ymm12
  953984:	vmulps ymm12,ymm4,ymm9
  953989:	vmulps ymm4,ymm4,ymm10
  95398e:	vsubps ymm12,ymm12,ymm14
  953993:	vmulps ymm14,ymm9,ymm15
  953998:	vaddps ymm4,ymm4,ymm14
  95399d:	vmulps ymm14,ymm3,ymm11
  9539a2:	vmulps ymm3,ymm8,ymm3
  9539a6:	vmovaps YMMWORD PTR [rsp+0xe0],ymm4
  9539af:	vaddps ymm4,ymm13,ymm12
  9539b4:	vmulps ymm12,ymm2,ymm8
  9539b9:	vmulps ymm2,ymm2,ymm11
  9539be:	vsubps ymm12,ymm14,ymm12
  9539c3:	vaddps ymm2,ymm3,ymm2
  9539c7:	vaddps ymm3,ymm13,ymm12
  9539cc:	cmp    rbx,r15
  9539cf:	jne    953880 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x260>
  9539d5:	lea    r15,[rcx+r10*1]
  9539d9:	vmovaps YMMWORD PTR [r15],ymm1
  9539de:	mov    rax,QWORD PTR [r14]
  9539e1:	add    r10,0x20
  9539e5:	vmovaps ymm6,YMMWORD PTR [rsp+0xa0]
  9539ee:	vmovaps ymm7,YMMWORD PTR [rsp+0x100]
  9539f7:	vmovaps YMMWORD PTR [rax+r11*1],ymm6
  9539fd:	mov    rax,QWORD PTR [r14]
  953a00:	add    r11,0x20
  953a04:	vmovaps ymm6,YMMWORD PTR [rsp+0x120]
  953a0d:	vmovaps YMMWORD PTR [rax+rdi*1],ymm0
  953a12:	mov    rax,QWORD PTR [r14]
  953a15:	mov    rdi,QWORD PTR [rsp+0x50]
  953a1a:	vmovaps YMMWORD PTR [rax+rsi*1],ymm7
  953a1f:	mov    rax,QWORD PTR [r14]
  953a22:	vmovaps ymm7,YMMWORD PTR [rsp+0xc0]
  953a2b:	mov    rsi,QWORD PTR [rsp+0x58]
  953a30:	vmovaps YMMWORD PTR [rax+r9*1],ymm6
  953a36:	mov    rax,QWORD PTR [r14]
  953a39:	vmovaps ymm6,YMMWORD PTR [rsp+0x80]
  953a42:	vmovaps YMMWORD PTR [rax+r8*1],ymm7
  953a48:	mov    rax,QWORD PTR [r14]
  953a4b:	vmovaps ymm7,YMMWORD PTR [rsp+0xe0]
  953a54:	vmovaps YMMWORD PTR [rax+r13*1],ymm5
  953a5a:	mov    rax,QWORD PTR [r14]
  953a5d:	vmovaps YMMWORD PTR [rax+r12*1],ymm6
  953a63:	mov    rax,QWORD PTR [r14]
  953a66:	vmovaps YMMWORD PTR [rax+rsi*1],ymm4
  953a6b:	mov    rax,QWORD PTR [r14]
  953a6e:	mov    rsi,QWORD PTR [rsp+0x48]
  953a73:	vmovaps YMMWORD PTR [rax+rdi*1],ymm7
  953a78:	mov    rax,QWORD PTR [r14]
  953a7b:	mov    rdi,QWORD PTR [rsp+0x40]
  953a80:	vmovaps YMMWORD PTR [rax+rsi*1],ymm3
  953a85:	mov    rax,QWORD PTR [r14]
  953a88:	vmovaps YMMWORD PTR [rax+rdi*1],ymm2
  953a8d:	mov    rsi,QWORD PTR [rsp+0x70]
  953a92:	add    QWORD PTR [rsp+0x78],0x8
  953a98:	mov    rax,QWORD PTR [rsp+0x78]
  953a9d:	cmp    rax,QWORD PTR [rsi+0x40]
  953aa1:	jb     953758 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x138>
  953aa7:	vzeroupper
  953aaa:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  953ab6:	lea    r12,[rax-0x360]
  953abd:	mov    r13,rax
  953ac0:	mov    r14,QWORD PTR [r12+0x188]
  953ac8:	test   r14,r14
  953acb:	je     953aef <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4cf>
  953acd:	movzx  eax,BYTE PTR [r12+0x180]
  953ad6:	lea    ebx,[rax-0x1]
  953ad9:	mov    BYTE PTR [r12+0x180],bl
  953ae1:	cmp    bl,0xf
  953ae4:	jbe    953afe <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4de>
  953ae6:	add    WORD PTR [r14+0xfe2],0x1
  953aef:	lea    rsp,[rbp-0x28]
  953af3:	pop    rbx
  953af4:	pop    r12
  953af6:	pop    r13
  953af8:	pop    r14
  953afa:	pop    r15
  953afc:	pop    rbp
  953afd:	ret
  953afe:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  953b03:	movzx  r15d,bl
  953b07:	movsxd rcx,r15d
  953b0a:	mov    rsi,rax
  953b0d:	vmovq  xmm6,rax
  953b12:	lea    rdx,[rcx+rcx*2]
  953b16:	shl    rdx,0x3
  953b1a:	sub    rsi,QWORD PTR [r12+rdx*1]
  953b1e:	vpinsrq xmm0,xmm6,rsi,0x1
  953b24:	vpsubq xmm0,xmm0,XMMWORD PTR [r13+rdx*1-0x360]
  953b2e:	movzx  r13d,BYTE PTR [r12+0x181]
  953b37:	cmp    r13b,bl
  953b3a:	jae    953b50 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x530>
  953b3c:	lea    eax,[r15-0x1]
  953b40:	vmovq  rdx,xmm0
  953b45:	cdqe
  953b47:	lea    rax,[rax+rax*2]
  953b4b:	add    QWORD PTR [r12+rax*8+0x8],rdx
  953b50:	movzx  r12d,WORD PTR [r14+0xfe0]
  953b58:	cmp    r12w,0x7f
  953b5d:	je     953cf0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x6d0>
  953b63:	mov    QWORD PTR [rsp+0x100],rcx
  953b6b:	sub    ebx,r13d
  953b6e:	vmovdqa XMMWORD PTR [rsp+0x120],xmm0
  953b77:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  953b83:	mov    rcx,QWORD PTR [rsp+0x100]
  953b8b:	mov    DWORD PTR [rsp+0x158],0x0
  953b96:	vmovdqa xmm0,XMMWORD PTR [rsp+0x120]
  953b9f:	vmovdqu XMMWORD PTR [rsp+0x148],xmm0
  953ba8:	lea    rdx,[rcx+rcx*2]
  953bac:	lea    rdx,[rax+rdx*8-0x360]
  953bb4:	mov    rax,QWORD PTR [rdx]
  953bb7:	mov    QWORD PTR [rsp+0x140],rax
  953bbf:	movzx  eax,BYTE PTR [rdx+0x10]
  953bc3:	lea    edx,[r12+0x1]
  953bc8:	vmovdqa xmm7,XMMWORD PTR [rsp+0x140]
  953bd1:	mov    WORD PTR [r14+0xfe0],dx
  953bd9:	mov    ah,bl
  953bdb:	mov    WORD PTR [rsp+0x15c],ax
  953be3:	movzx  eax,r12w
  953be7:	vmovdqu xmm6,XMMWORD PTR [rsp+0x14e]
  953bf0:	shl    rax,0x5
  953bf4:	add    rax,r14
  953bf7:	vmovdqu XMMWORD PTR [rax],xmm7
  953bfb:	vmovdqu XMMWORD PTR [rax+0xe],xmm6
  953c00:	lea    rsp,[rbp-0x28]
  953c04:	pop    rbx
  953c05:	pop    r12
  953c07:	pop    r13
  953c09:	pop    r14
  953c0b:	pop    r15
  953c0d:	pop    rbp
  953c0e:	ret
  953c0f:	shl    r8,0x4
  953c13:	add    rdi,r8
  953c16:	cmp    r9,rdi
  953c19:	je     9536e2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc2>
  953c1f:	mov    QWORD PTR [r14+0x8],rdi
  953c23:	jmp    9536e2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc2>
  953c28:	mov    eax,0x1
  953c2d:	lock xadd DWORD PTR [rip+0x2d2f2b],eax        # c26b60 <cpl::Profiling::registerRegion(char const*)::counter>
  953c35:	mov    edx,0x1
  953c3a:	mov    r13d,0x1
  953c40:	add    eax,0x2
  953c43:	cmp    eax,0xfe
  953c48:	ja     953c63 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x643>
  953c4a:	lea    rdx,[rip+0x2d13cf]        # c25020 <cpl::Profiling::regions>
  953c51:	mov    ecx,eax
  953c53:	mov    r13d,eax
  953c56:	lea    rsi,[rip+0xac655]        # a002b2 <_IO_stdin_used+0x182b2>
  953c5d:	mov    QWORD PTR [rdx+rcx*8],rsi
  953c61:	mov    edx,eax
  953c63:	xor    eax,eax
  953c65:	lock cmpxchg BYTE PTR [rip+0x2de491],dl        # c320fe <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  953c6d:	cmovne r13d,eax
  953c71:	jmp    95365d <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x3d>
  953c76:	sub    rsi,rcx
  953c79:	mov    rdi,r14
  953c7c:	call   939b50 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  953c81:	mov    rdi,QWORD PTR [rsp+0x70]
  953c86:	mov    rax,QWORD PTR [rdi+0x38]
  953c8a:	mov    rdx,QWORD PTR [rdi+0x48]
  953c8e:	jmp    9536e2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc2>
  953c93:	mov    QWORD PTR [rsp+0x120],rax
  953c9b:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  953ca0:	mov    rdx,QWORD PTR [rsp+0x120]
  953ca8:	mov    QWORD PTR [rsp+0x148],0x0
  953cb4:	mov    QWORD PTR [rsp+0x140],rax
  953cbc:	lea    rax,[r12+r12*2]
  953cc0:	vmovdqa xmm7,XMMWORD PTR [rsp+0x140]
  953cc9:	shl    rax,0x3
  953ccd:	mov    BYTE PTR [rsp+0x150],r13b
  953cd5:	vmovdqu XMMWORD PTR [rdx+rax*1-0x360],xmm7
  953cde:	mov    BYTE PTR [rax+r15*1+0x10],r13b
  953ce3:	movzx  r12d,BYTE PTR [r15+0x180]
  953ceb:	jmp    95368c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x6c>
  953cf0:	movzx  eax,WORD PTR [r14+0xfe2]
  953cf8:	cmp    ax,0xffff
  953cfc:	je     953aef <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4cf>
  953d02:	add    eax,0x1
  953d05:	mov    WORD PTR [r14+0xfe2],ax
  953d0d:	jmp    953aef <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(8), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4cf>
  953d12:	endbr64
  953d16:	mov    rbx,rax
  953d19:	vzeroupper
  953d1c:	call   93c730 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  953d21:	mov    rdi,rbx
  953d24:	call   e1b10 <_Unwind_Resume@plt>

Disassembly of section .fini:
