; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000958860 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long):
  958860:	endbr64
  958864:	push   r15
  958866:	mov    r15,rcx
  958869:	push   r14
  95886b:	push   r13
  95886d:	push   r12
  95886f:	push   rbp
  958870:	push   rbx
  958871:	mov    rbx,rdi
  958874:	sub    rsp,0xf8
  95887b:	mov    QWORD PTR [rsp+0x70],rsi
  958880:	mov    QWORD PTR [rsp+0xa0],rdx
  958888:	movzx  r13d,BYTE PTR [rip+0x2d9862]        # c320f2 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  958890:	test   r13b,r13b
  958893:	je     958e05 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x5a5>
  958899:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  9588a5:	cmp    QWORD PTR [rax-0x1d8],0x0
  9588ad:	mov    rbp,rax
  9588b0:	lea    r14,[rax-0x360]
  9588b7:	je     9588e2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x82>
  9588b9:	movzx  r12d,BYTE PTR [r14+0x180]
  9588c1:	cmp    r12b,0xf
  9588c5:	jbe    958e70 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x610>
  9588cb:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  9588d7:	add    r12d,0x1
  9588db:	mov    BYTE PTR [rax-0x1e0],r12b
  9588e2:	mov    rax,QWORD PTR [rsp+0x70]
  9588e7:	mov    r9,QWORD PTR [rbx+0x8]
  9588eb:	mov    rdi,QWORD PTR [rbx]
  9588ee:	mov    rdx,QWORD PTR [rax+0x48]
  9588f2:	mov    rax,QWORD PTR [rax+0x38]
  9588f6:	mov    rcx,r9
  9588f9:	sub    rcx,rdi
  9588fc:	sar    rcx,0x2
  958900:	mov    r8,rax
  958903:	imul   r8,rdx
  958907:	lea    rsi,[r8*4+0x0]
  95890f:	cmp    rcx,rsi
  958912:	jb     958e53 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x5f3>
  958918:	cmp    rsi,rcx
  95891b:	jb     958dec <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x58c>
  958921:	mov    rsi,QWORD PTR [rsp+0x70]
  958926:	lea    rcx,[rdx+rdx*1]
  95892a:	imul   rax,rcx
  95892e:	cmp    QWORD PTR [rsi+0x40],0x0
  958933:	je     958cb2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x452>
  958939:	lea    r11,[rdx*4+0x0]
  958941:	xor    r10d,r10d
  958944:	mov    QWORD PTR [rsp+0x68],0x0
  95894d:	lea    rsi,[r11+rax*1]
  958951:	shl    rax,0x2
  958955:	mov    QWORD PTR [rsp+0xa8],r11
  95895d:	mov    QWORD PTR [rsp+0xc0],rax
  958965:	mov    rax,rsi
  958968:	lea    rdi,[rdx*8+0x0]
  958970:	sub    rax,rcx
  958973:	shl    rdx,0x4
  958977:	shl    rax,0x2
  95897b:	mov    QWORD PTR [rsp+0xb0],rdi
  958983:	mov    QWORD PTR [rsp+0xb8],rax
  95898b:	lea    rax,[rsi*4+0x0]
  958993:	mov    QWORD PTR [rsp+0xc8],rdx
  95899b:	mov    QWORD PTR [rsp+0x78],rax
  9589a0:	mov    rdi,QWORD PTR [rsp+0xb0]
  9589a8:	mov    rcx,QWORD PTR [rsp+0xa8]
  9589b0:	mov    rax,QWORD PTR [rsp+0x70]
  9589b5:	mov    rbp,QWORD PTR [rsp+0xc0]
  9589bd:	add    rdi,r10
  9589c0:	mov    rdx,QWORD PTR [rsp+0xb8]
  9589c8:	lea    rsi,[rcx+rdi*1]
  9589cc:	lea    r9,[rcx+rsi*1]
  9589d0:	mov    rcx,QWORD PTR [rsp+0xc8]
  9589d8:	mov    rax,QWORD PTR [rax]
  9589db:	lea    r13,[rbp+r10*1+0x0]
  9589e0:	lea    r12,[rbp+r11*1+0x0]
  9589e5:	mov    rbp,rdx
  9589e8:	add    rdx,r11
  9589eb:	add    rbp,r10
  9589ee:	mov    QWORD PTR [rsp+0x80],rdx
  9589f6:	lea    r8,[rcx+r11*1]
  9589fa:	mov    rcx,QWORD PTR [rsp+0xa0]
  958a02:	vmovaps xmm6,XMMWORD PTR [rax+r10*1]
  958a08:	vmovaps xmm7,XMMWORD PTR [rax+r11*1]
  958a0e:	vmovaps xmm9,XMMWORD PTR [rax+rdi*1]
  958a13:	vmovaps xmm10,XMMWORD PTR [rax+rsi*1]
  958a18:	vmovaps xmm11,XMMWORD PTR [rax+r9*1]
  958a1e:	vmovaps xmm8,XMMWORD PTR [rax+r8*1]
  958a24:	mov    rax,QWORD PTR [rcx]
  958a27:	mov    rcx,QWORD PTR [rcx+0x10]
  958a2b:	mov    QWORD PTR [rsp+0x98],rax
  958a33:	mov    rax,QWORD PTR [rbx]
  958a36:	vmovaps xmm0,XMMWORD PTR [rax+r8*1]
  958a3c:	lea    r14,[rax+r10*1]
  958a40:	vmovaps xmm5,XMMWORD PTR [r14]
  958a45:	vmovaps XMMWORD PTR [rsp+0x30],xmm0
  958a4b:	vmovaps xmm0,XMMWORD PTR [rax+r13*1]
  958a51:	vmovaps xmm4,XMMWORD PTR [rax+r11*1]
  958a57:	vmovaps XMMWORD PTR [rsp+0x10],xmm0
  958a5d:	vmovaps xmm0,XMMWORD PTR [rax+r12*1]
  958a63:	vmovaps xmm3,XMMWORD PTR [rax+rdi*1]
  958a68:	vmovaps XMMWORD PTR [rsp+0x50],xmm0
  958a6e:	vmovaps xmm0,XMMWORD PTR [rax+rbp*1]
  958a73:	vmovaps xmm2,XMMWORD PTR [rax+rsi*1]
  958a78:	vmovaps xmm1,XMMWORD PTR [rax+r9*1]
  958a7e:	vmovaps XMMWORD PTR [rsp],xmm0
  958a83:	vmovaps xmm0,XMMWORD PTR [rax+rdx*1]
  958a88:	mov    rdx,QWORD PTR [rsp+0x78]
  958a8d:	vmovaps XMMWORD PTR [rsp+0x40],xmm0
  958a93:	add    rdx,r10
  958a96:	vmovaps xmm0,XMMWORD PTR [rax+rdx*1]
  958a9b:	mov    QWORD PTR [rsp+0x88],rdx
  958aa3:	mov    rdx,QWORD PTR [rsp+0x78]
  958aa8:	vmovaps XMMWORD PTR [rsp+0x20],xmm0
  958aae:	add    rdx,r11
  958ab1:	vmovaps xmm0,XMMWORD PTR [rax+rdx*1]
  958ab6:	mov    QWORD PTR [rsp+0x90],rdx
  958abe:	test   r15,r15
  958ac1:	je     958bf3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x393>
  958ac7:	mov    rdx,QWORD PTR [rsp+0x98]
  958acf:	xor    eax,eax
  958ad1:	nop    DWORD PTR [rax+0x0]
  958ad8:	vmulps xmm12,xmm7,xmm4
  958adc:	vbroadcastss xmm13,DWORD PTR [rdx+rax*4]
  958ae2:	vmulps xmm14,xmm5,xmm7
  958ae6:	vmovaps xmm15,XMMWORD PTR [rsp+0x30]
  958aec:	vmulps xmm5,xmm5,xmm6
  958af0:	vmulps xmm4,xmm6,xmm4
  958af4:	vsubps xmm12,xmm13,xmm12
  958af9:	vaddps xmm5,xmm5,xmm12
  958afe:	vmulps xmm12,xmm2,xmm10
  958b03:	vaddps xmm4,xmm14,xmm4
  958b07:	vmulps xmm14,xmm10,xmm3
  958b0b:	vmulps xmm3,xmm9,xmm3
  958b0f:	vmulps xmm2,xmm2,xmm9
  958b14:	vsubps xmm12,xmm13,xmm12
  958b19:	vaddps xmm3,xmm3,xmm12
  958b1e:	vmulps xmm12,xmm8,xmm15
  958b23:	vaddps xmm2,xmm14,xmm2
  958b27:	vmulps xmm14,xmm8,xmm1
  958b2b:	vmulps xmm1,xmm11,xmm1
  958b2f:	vsubps xmm12,xmm13,xmm12
  958b34:	vmulps xmm13,xmm11,xmm15
  958b39:	vmulps xmm15,xmm7,XMMWORD PTR [rsp+0x10]
  958b3f:	vaddps xmm1,xmm1,xmm12
  958b44:	vbroadcastss xmm12,DWORD PTR [rcx+rax*4]
  958b4a:	add    rax,0x1
  958b4e:	vaddps xmm13,xmm14,xmm13
  958b53:	vmovaps xmm14,XMMWORD PTR [rsp+0x50]
  958b59:	vmovaps XMMWORD PTR [rsp+0x30],xmm13
  958b5f:	vmulps xmm13,xmm7,xmm14
  958b64:	vmulps xmm14,xmm6,xmm14
  958b69:	vsubps xmm13,xmm12,xmm13
  958b6e:	vaddps xmm14,xmm15,xmm14
  958b73:	vmulps xmm15,xmm10,XMMWORD PTR [rsp]
  958b78:	vmovaps XMMWORD PTR [rsp+0x50],xmm14
  958b7e:	vmulps xmm14,xmm6,XMMWORD PTR [rsp+0x10]
  958b84:	vaddps xmm13,xmm14,xmm13
  958b89:	vmovaps xmm14,XMMWORD PTR [rsp+0x40]
  958b8f:	vmovaps XMMWORD PTR [rsp+0x10],xmm13
  958b95:	vmulps xmm13,xmm10,xmm14
  958b9a:	vmulps xmm14,xmm9,xmm14
  958b9f:	vsubps xmm13,xmm12,xmm13
  958ba4:	vaddps xmm14,xmm15,xmm14
  958ba9:	vmovaps xmm15,XMMWORD PTR [rsp+0x20]
  958baf:	vmovaps XMMWORD PTR [rsp+0x40],xmm14
  958bb5:	vmulps xmm14,xmm9,XMMWORD PTR [rsp]
  958bba:	vaddps xmm13,xmm14,xmm13
  958bbf:	vmovaps XMMWORD PTR [rsp],xmm13
  958bc4:	vmulps xmm13,xmm8,xmm0
  958bc8:	vmulps xmm0,xmm11,xmm0
  958bcc:	vsubps xmm12,xmm12,xmm13
  958bd1:	vmulps xmm13,xmm8,xmm15
  958bd6:	vaddps xmm0,xmm13,xmm0
  958bda:	vmulps xmm13,xmm11,xmm15
  958bdf:	vaddps xmm12,xmm13,xmm12
  958be4:	vmovaps XMMWORD PTR [rsp+0x20],xmm12
  958bea:	cmp    r15,rax
  958bed:	jne    958ad8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x278>
  958bf3:	vmovaps XMMWORD PTR [r14],xmm5
  958bf8:	mov    rax,QWORD PTR [rbx]
  958bfb:	add    r10,0x10
  958bff:	vmovaps xmm6,XMMWORD PTR [rsp+0x30]
  958c05:	vmovaps xmm7,XMMWORD PTR [rsp+0x10]
  958c0b:	vmovaps XMMWORD PTR [rax+r11*1],xmm4
  958c11:	mov    rax,QWORD PTR [rbx]
  958c14:	add    r11,0x10
  958c18:	vmovaps XMMWORD PTR [rax+rdi*1],xmm3
  958c1d:	mov    rax,QWORD PTR [rbx]
  958c20:	mov    rdi,QWORD PTR [rsp+0x88]
  958c28:	vmovaps XMMWORD PTR [rax+rsi*1],xmm2
  958c2d:	mov    rax,QWORD PTR [rbx]
  958c30:	mov    rsi,QWORD PTR [rsp+0x80]
  958c38:	vmovaps XMMWORD PTR [rax+r9*1],xmm1
  958c3e:	mov    rax,QWORD PTR [rbx]
  958c41:	vmovaps XMMWORD PTR [rax+r8*1],xmm6
  958c47:	mov    rax,QWORD PTR [rbx]
  958c4a:	vmovaps xmm6,XMMWORD PTR [rsp+0x50]
  958c50:	vmovaps XMMWORD PTR [rax+r13*1],xmm7
  958c56:	mov    rax,QWORD PTR [rbx]
  958c59:	vmovaps xmm7,XMMWORD PTR [rsp]
  958c5e:	vmovaps XMMWORD PTR [rax+r12*1],xmm6
  958c64:	mov    rax,QWORD PTR [rbx]
  958c67:	vmovaps xmm6,XMMWORD PTR [rsp+0x40]
  958c6d:	vmovaps XMMWORD PTR [rax+rbp*1],xmm7
  958c72:	mov    rax,QWORD PTR [rbx]
  958c75:	vmovaps xmm7,XMMWORD PTR [rsp+0x20]
  958c7b:	vmovaps XMMWORD PTR [rax+rsi*1],xmm6
  958c80:	mov    rax,QWORD PTR [rbx]
  958c83:	mov    rsi,QWORD PTR [rsp+0x90]
  958c8b:	vmovaps XMMWORD PTR [rax+rdi*1],xmm7
  958c90:	mov    rax,QWORD PTR [rbx]
  958c93:	vmovaps XMMWORD PTR [rax+rsi*1],xmm0
  958c98:	add    QWORD PTR [rsp+0x68],0x4
  958c9e:	mov    rdi,QWORD PTR [rsp+0x70]
  958ca3:	mov    rax,QWORD PTR [rsp+0x68]
  958ca8:	cmp    rax,QWORD PTR [rdi+0x40]
  958cac:	jb     9589a0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x140>
  958cb2:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  958cbe:	lea    rbp,[rax-0x360]
  958cc5:	mov    r12,rax
  958cc8:	mov    r13,QWORD PTR [rbp+0x188]
  958ccf:	test   r13,r13
  958cd2:	je     958cf2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x492>
  958cd4:	movzx  eax,BYTE PTR [rbp+0x180]
  958cdb:	lea    ebx,[rax-0x1]
  958cde:	mov    BYTE PTR [rbp+0x180],bl
  958ce4:	cmp    bl,0xf
  958ce7:	jbe    958d04 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4a4>
  958ce9:	add    WORD PTR [r13+0xfe2],0x1
  958cf2:	add    rsp,0xf8
  958cf9:	pop    rbx
  958cfa:	pop    rbp
  958cfb:	pop    r12
  958cfd:	pop    r13
  958cff:	pop    r14
  958d01:	pop    r15
  958d03:	ret
  958d04:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  958d09:	movzx  r15d,bl
  958d0d:	lea    rdx,[r15+r15*2]
  958d11:	mov    rcx,rax
  958d14:	vmovq  xmm6,rax
  958d19:	shl    rdx,0x3
  958d1d:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  958d22:	vpinsrq xmm0,xmm6,rcx,0x1
  958d28:	vpsubq xmm0,xmm0,XMMWORD PTR [r12+rdx*1-0x360]
  958d32:	movzx  r12d,BYTE PTR [rbp+0x181]
  958d3a:	cmp    r12b,bl
  958d3d:	jae    958d53 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4f3>
  958d3f:	lea    eax,[r15-0x1]
  958d43:	vmovq  rdx,xmm0
  958d48:	cdqe
  958d4a:	lea    rax,[rax+rax*2]
  958d4e:	add    QWORD PTR [rbp+rax*8+0x8],rdx
  958d53:	movzx  ebp,WORD PTR [r13+0xfe0]
  958d5b:	cmp    bp,0x7f
  958d5f:	je     958ebd <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x65d>
  958d65:	vmovdqa XMMWORD PTR [rsp],xmm0
  958d6a:	sub    ebx,r12d
  958d6d:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  958d79:	lea    rdx,[r15+r15*2]
  958d7d:	vmovdqa xmm0,XMMWORD PTR [rsp]
  958d82:	mov    DWORD PTR [rsp+0xe8],0x0
  958d8d:	vmovdqu XMMWORD PTR [rsp+0xd8],xmm0
  958d96:	lea    rdx,[rax+rdx*8-0x360]
  958d9e:	mov    rax,QWORD PTR [rdx]
  958da1:	mov    QWORD PTR [rsp+0xd0],rax
  958da9:	movzx  eax,BYTE PTR [rdx+0x10]
  958dad:	lea    edx,[rbp+0x1]
  958db0:	vmovdqa xmm7,XMMWORD PTR [rsp+0xd0]
  958db9:	mov    WORD PTR [r13+0xfe0],dx
  958dc1:	mov    ah,bl
  958dc3:	mov    WORD PTR [rsp+0xec],ax
  958dcb:	movzx  eax,bp
  958dce:	vmovdqu xmm6,XMMWORD PTR [rsp+0xde]
  958dd7:	shl    rax,0x5
  958ddb:	add    rax,r13
  958dde:	vmovdqu XMMWORD PTR [rax],xmm7
  958de2:	vmovdqu XMMWORD PTR [rax+0xe],xmm6
  958de7:	jmp    958cf2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x492>
  958dec:	shl    r8,0x4
  958df0:	add    rdi,r8
  958df3:	cmp    r9,rdi
  958df6:	je     958921 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc1>
  958dfc:	mov    QWORD PTR [rbx+0x8],rdi
  958e00:	jmp    958921 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc1>
  958e05:	mov    eax,0x1
  958e0a:	lock xadd DWORD PTR [rip+0x2cdd4e],eax        # c26b60 <cpl::Profiling::registerRegion(char const*)::counter>
  958e12:	mov    edx,0x1
  958e17:	mov    r13d,0x1
  958e1d:	add    eax,0x2
  958e20:	cmp    eax,0xfe
  958e25:	ja     958e40 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x5e0>
  958e27:	lea    rdx,[rip+0x2cc1f2]        # c25020 <cpl::Profiling::regions>
  958e2e:	mov    ecx,eax
  958e30:	mov    r13d,eax
  958e33:	lea    rsi,[rip+0xa7478]        # a002b2 <_IO_stdin_used+0x182b2>
  958e3a:	mov    QWORD PTR [rdx+rcx*8],rsi
  958e3e:	mov    edx,eax
  958e40:	xor    eax,eax
  958e42:	lock cmpxchg BYTE PTR [rip+0x2d92a8],dl        # c320f2 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  958e4a:	cmovne r13d,eax
  958e4e:	jmp    958899 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x39>
  958e53:	sub    rsi,rcx
  958e56:	mov    rdi,rbx
  958e59:	call   939b50 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  958e5e:	mov    rsi,QWORD PTR [rsp+0x70]
  958e63:	mov    rax,QWORD PTR [rsi+0x38]
  958e67:	mov    rdx,QWORD PTR [rsi+0x48]
  958e6b:	jmp    958921 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc1>
  958e70:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  958e75:	mov    BYTE PTR [rsp+0xe0],r13b
  958e7d:	mov    QWORD PTR [rsp+0xd0],rax
  958e85:	lea    rax,[r12+r12*2]
  958e89:	mov    QWORD PTR [rsp+0xd8],0x0
  958e95:	shl    rax,0x3
  958e99:	vmovdqa xmm7,XMMWORD PTR [rsp+0xd0]
  958ea2:	vmovdqu XMMWORD PTR [rbp+rax*1-0x360],xmm7
  958eab:	mov    BYTE PTR [rax+r14*1+0x10],r13b
  958eb0:	movzx  r12d,BYTE PTR [r14+0x180]
  958eb8:	jmp    9588cb <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x6b>
  958ebd:	movzx  eax,WORD PTR [r13+0xfe2]
  958ec5:	cmp    ax,0xffff
  958ec9:	je     958cf2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x492>
  958ecf:	add    eax,0x1
  958ed2:	mov    WORD PTR [r13+0xfe2],ax
  958eda:	jmp    958cf2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x492>
  958edf:	endbr64
  958ee3:	mov    rbx,rax
  958ee6:	vzeroupper
  958ee9:	call   93c730 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  958eee:	mov    rdi,rbx
  958ef1:	call   e1b10 <_Unwind_Resume@plt>

Disassembly of section .fini:
