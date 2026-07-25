; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

000000000095ca40 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long):
  95ca40:	endbr64
  95ca44:	push   r15
  95ca46:	mov    r15,rcx
  95ca49:	push   r14
  95ca4b:	push   r13
  95ca4d:	push   r12
  95ca4f:	push   rbp
  95ca50:	push   rbx
  95ca51:	mov    rbx,rdi
  95ca54:	sub    rsp,0xf8
  95ca5b:	mov    QWORD PTR [rsp+0x70],rsi
  95ca60:	mov    QWORD PTR [rsp+0xa0],rdx
  95ca68:	movzx  r13d,BYTE PTR [rip+0x2d5678]        # c320e8 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  95ca70:	test   r13b,r13b
  95ca73:	je     95cfe5 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x5a5>
  95ca79:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  95ca85:	cmp    QWORD PTR [rax-0x1d8],0x0
  95ca8d:	mov    rbp,rax
  95ca90:	lea    r14,[rax-0x360]
  95ca97:	je     95cac2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x82>
  95ca99:	movzx  r12d,BYTE PTR [r14+0x180]
  95caa1:	cmp    r12b,0xf
  95caa5:	jbe    95d050 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x610>
  95caab:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  95cab7:	add    r12d,0x1
  95cabb:	mov    BYTE PTR [rax-0x1e0],r12b
  95cac2:	mov    rax,QWORD PTR [rsp+0x70]
  95cac7:	mov    r9,QWORD PTR [rbx+0x8]
  95cacb:	mov    rdi,QWORD PTR [rbx]
  95cace:	mov    rdx,QWORD PTR [rax+0x48]
  95cad2:	mov    rax,QWORD PTR [rax+0x38]
  95cad6:	mov    rcx,r9
  95cad9:	sub    rcx,rdi
  95cadc:	sar    rcx,0x2
  95cae0:	mov    r8,rax
  95cae3:	imul   r8,rdx
  95cae7:	lea    rsi,[r8*4+0x0]
  95caef:	cmp    rcx,rsi
  95caf2:	jb     95d033 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x5f3>
  95caf8:	cmp    rsi,rcx
  95cafb:	jb     95cfcc <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x58c>
  95cb01:	mov    rsi,QWORD PTR [rsp+0x70]
  95cb06:	lea    rcx,[rdx+rdx*1]
  95cb0a:	imul   rax,rcx
  95cb0e:	cmp    QWORD PTR [rsi+0x40],0x0
  95cb13:	je     95ce92 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x452>
  95cb19:	lea    r11,[rdx*4+0x0]
  95cb21:	xor    r10d,r10d
  95cb24:	mov    QWORD PTR [rsp+0x68],0x0
  95cb2d:	lea    rsi,[r11+rax*1]
  95cb31:	shl    rax,0x2
  95cb35:	mov    QWORD PTR [rsp+0xa8],r11
  95cb3d:	mov    QWORD PTR [rsp+0xc0],rax
  95cb45:	mov    rax,rsi
  95cb48:	lea    rdi,[rdx*8+0x0]
  95cb50:	sub    rax,rcx
  95cb53:	shl    rdx,0x4
  95cb57:	shl    rax,0x2
  95cb5b:	mov    QWORD PTR [rsp+0xb0],rdi
  95cb63:	mov    QWORD PTR [rsp+0xb8],rax
  95cb6b:	lea    rax,[rsi*4+0x0]
  95cb73:	mov    QWORD PTR [rsp+0xc8],rdx
  95cb7b:	mov    QWORD PTR [rsp+0x78],rax
  95cb80:	mov    rdi,QWORD PTR [rsp+0xb0]
  95cb88:	mov    rcx,QWORD PTR [rsp+0xa8]
  95cb90:	mov    rax,QWORD PTR [rsp+0x70]
  95cb95:	mov    rbp,QWORD PTR [rsp+0xc0]
  95cb9d:	add    rdi,r10
  95cba0:	mov    rdx,QWORD PTR [rsp+0xb8]
  95cba8:	lea    rsi,[rcx+rdi*1]
  95cbac:	lea    r9,[rcx+rsi*1]
  95cbb0:	mov    rcx,QWORD PTR [rsp+0xc8]
  95cbb8:	mov    rax,QWORD PTR [rax]
  95cbbb:	lea    r13,[rbp+r10*1+0x0]
  95cbc0:	lea    r12,[rbp+r11*1+0x0]
  95cbc5:	mov    rbp,rdx
  95cbc8:	add    rdx,r11
  95cbcb:	add    rbp,r10
  95cbce:	mov    QWORD PTR [rsp+0x80],rdx
  95cbd6:	lea    r8,[rcx+r11*1]
  95cbda:	mov    rcx,QWORD PTR [rsp+0xa0]
  95cbe2:	vmovaps xmm6,XMMWORD PTR [rax+r10*1]
  95cbe8:	vmovaps xmm7,XMMWORD PTR [rax+r11*1]
  95cbee:	vmovaps xmm9,XMMWORD PTR [rax+rdi*1]
  95cbf3:	vmovaps xmm10,XMMWORD PTR [rax+rsi*1]
  95cbf8:	vmovaps xmm11,XMMWORD PTR [rax+r9*1]
  95cbfe:	vmovaps xmm8,XMMWORD PTR [rax+r8*1]
  95cc04:	mov    rax,QWORD PTR [rcx]
  95cc07:	mov    rcx,QWORD PTR [rcx+0x10]
  95cc0b:	mov    QWORD PTR [rsp+0x98],rax
  95cc13:	mov    rax,QWORD PTR [rbx]
  95cc16:	vmovaps xmm0,XMMWORD PTR [rax+r8*1]
  95cc1c:	lea    r14,[rax+r10*1]
  95cc20:	vmovaps xmm5,XMMWORD PTR [r14]
  95cc25:	vmovaps XMMWORD PTR [rsp+0x30],xmm0
  95cc2b:	vmovaps xmm0,XMMWORD PTR [rax+r13*1]
  95cc31:	vmovaps xmm4,XMMWORD PTR [rax+r11*1]
  95cc37:	vmovaps XMMWORD PTR [rsp+0x10],xmm0
  95cc3d:	vmovaps xmm0,XMMWORD PTR [rax+r12*1]
  95cc43:	vmovaps xmm3,XMMWORD PTR [rax+rdi*1]
  95cc48:	vmovaps XMMWORD PTR [rsp+0x50],xmm0
  95cc4e:	vmovaps xmm0,XMMWORD PTR [rax+rbp*1]
  95cc53:	vmovaps xmm2,XMMWORD PTR [rax+rsi*1]
  95cc58:	vmovaps xmm1,XMMWORD PTR [rax+r9*1]
  95cc5e:	vmovaps XMMWORD PTR [rsp],xmm0
  95cc63:	vmovaps xmm0,XMMWORD PTR [rax+rdx*1]
  95cc68:	mov    rdx,QWORD PTR [rsp+0x78]
  95cc6d:	vmovaps XMMWORD PTR [rsp+0x40],xmm0
  95cc73:	add    rdx,r10
  95cc76:	vmovaps xmm0,XMMWORD PTR [rax+rdx*1]
  95cc7b:	mov    QWORD PTR [rsp+0x88],rdx
  95cc83:	mov    rdx,QWORD PTR [rsp+0x78]
  95cc88:	vmovaps XMMWORD PTR [rsp+0x20],xmm0
  95cc8e:	add    rdx,r11
  95cc91:	vmovaps xmm0,XMMWORD PTR [rax+rdx*1]
  95cc96:	mov    QWORD PTR [rsp+0x90],rdx
  95cc9e:	test   r15,r15
  95cca1:	je     95cdd3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x393>
  95cca7:	mov    rdx,QWORD PTR [rsp+0x98]
  95ccaf:	xor    eax,eax
  95ccb1:	nop    DWORD PTR [rax+0x0]
  95ccb8:	vmulps xmm12,xmm7,xmm4
  95ccbc:	vbroadcastss xmm13,DWORD PTR [rdx+rax*4]
  95ccc2:	vmulps xmm14,xmm5,xmm7
  95ccc6:	vmovaps xmm15,XMMWORD PTR [rsp+0x30]
  95cccc:	vmulps xmm5,xmm5,xmm6
  95ccd0:	vmulps xmm4,xmm6,xmm4
  95ccd4:	vsubps xmm12,xmm13,xmm12
  95ccd9:	vaddps xmm5,xmm5,xmm12
  95ccde:	vmulps xmm12,xmm2,xmm10
  95cce3:	vaddps xmm4,xmm14,xmm4
  95cce7:	vmulps xmm14,xmm10,xmm3
  95cceb:	vmulps xmm3,xmm9,xmm3
  95ccef:	vmulps xmm2,xmm2,xmm9
  95ccf4:	vsubps xmm12,xmm13,xmm12
  95ccf9:	vaddps xmm3,xmm3,xmm12
  95ccfe:	vmulps xmm12,xmm8,xmm15
  95cd03:	vaddps xmm2,xmm14,xmm2
  95cd07:	vmulps xmm14,xmm8,xmm1
  95cd0b:	vmulps xmm1,xmm11,xmm1
  95cd0f:	vsubps xmm12,xmm13,xmm12
  95cd14:	vmulps xmm13,xmm11,xmm15
  95cd19:	vmulps xmm15,xmm7,XMMWORD PTR [rsp+0x10]
  95cd1f:	vaddps xmm1,xmm1,xmm12
  95cd24:	vbroadcastss xmm12,DWORD PTR [rcx+rax*4]
  95cd2a:	add    rax,0x1
  95cd2e:	vaddps xmm13,xmm14,xmm13
  95cd33:	vmovaps xmm14,XMMWORD PTR [rsp+0x50]
  95cd39:	vmovaps XMMWORD PTR [rsp+0x30],xmm13
  95cd3f:	vmulps xmm13,xmm7,xmm14
  95cd44:	vmulps xmm14,xmm6,xmm14
  95cd49:	vsubps xmm13,xmm12,xmm13
  95cd4e:	vaddps xmm14,xmm15,xmm14
  95cd53:	vmulps xmm15,xmm10,XMMWORD PTR [rsp]
  95cd58:	vmovaps XMMWORD PTR [rsp+0x50],xmm14
  95cd5e:	vmulps xmm14,xmm6,XMMWORD PTR [rsp+0x10]
  95cd64:	vaddps xmm13,xmm14,xmm13
  95cd69:	vmovaps xmm14,XMMWORD PTR [rsp+0x40]
  95cd6f:	vmovaps XMMWORD PTR [rsp+0x10],xmm13
  95cd75:	vmulps xmm13,xmm10,xmm14
  95cd7a:	vmulps xmm14,xmm9,xmm14
  95cd7f:	vsubps xmm13,xmm12,xmm13
  95cd84:	vaddps xmm14,xmm15,xmm14
  95cd89:	vmovaps xmm15,XMMWORD PTR [rsp+0x20]
  95cd8f:	vmovaps XMMWORD PTR [rsp+0x40],xmm14
  95cd95:	vmulps xmm14,xmm9,XMMWORD PTR [rsp]
  95cd9a:	vaddps xmm13,xmm14,xmm13
  95cd9f:	vmovaps XMMWORD PTR [rsp],xmm13
  95cda4:	vmulps xmm13,xmm8,xmm0
  95cda8:	vmulps xmm0,xmm11,xmm0
  95cdac:	vsubps xmm12,xmm12,xmm13
  95cdb1:	vmulps xmm13,xmm8,xmm15
  95cdb6:	vaddps xmm0,xmm13,xmm0
  95cdba:	vmulps xmm13,xmm11,xmm15
  95cdbf:	vaddps xmm12,xmm13,xmm12
  95cdc4:	vmovaps XMMWORD PTR [rsp+0x20],xmm12
  95cdca:	cmp    r15,rax
  95cdcd:	jne    95ccb8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x278>
  95cdd3:	vmovaps XMMWORD PTR [r14],xmm5
  95cdd8:	mov    rax,QWORD PTR [rbx]
  95cddb:	add    r10,0x10
  95cddf:	vmovaps xmm6,XMMWORD PTR [rsp+0x30]
  95cde5:	vmovaps xmm7,XMMWORD PTR [rsp+0x10]
  95cdeb:	vmovaps XMMWORD PTR [rax+r11*1],xmm4
  95cdf1:	mov    rax,QWORD PTR [rbx]
  95cdf4:	add    r11,0x10
  95cdf8:	vmovaps XMMWORD PTR [rax+rdi*1],xmm3
  95cdfd:	mov    rax,QWORD PTR [rbx]
  95ce00:	mov    rdi,QWORD PTR [rsp+0x88]
  95ce08:	vmovaps XMMWORD PTR [rax+rsi*1],xmm2
  95ce0d:	mov    rax,QWORD PTR [rbx]
  95ce10:	mov    rsi,QWORD PTR [rsp+0x80]
  95ce18:	vmovaps XMMWORD PTR [rax+r9*1],xmm1
  95ce1e:	mov    rax,QWORD PTR [rbx]
  95ce21:	vmovaps XMMWORD PTR [rax+r8*1],xmm6
  95ce27:	mov    rax,QWORD PTR [rbx]
  95ce2a:	vmovaps xmm6,XMMWORD PTR [rsp+0x50]
  95ce30:	vmovaps XMMWORD PTR [rax+r13*1],xmm7
  95ce36:	mov    rax,QWORD PTR [rbx]
  95ce39:	vmovaps xmm7,XMMWORD PTR [rsp]
  95ce3e:	vmovaps XMMWORD PTR [rax+r12*1],xmm6
  95ce44:	mov    rax,QWORD PTR [rbx]
  95ce47:	vmovaps xmm6,XMMWORD PTR [rsp+0x40]
  95ce4d:	vmovaps XMMWORD PTR [rax+rbp*1],xmm7
  95ce52:	mov    rax,QWORD PTR [rbx]
  95ce55:	vmovaps xmm7,XMMWORD PTR [rsp+0x20]
  95ce5b:	vmovaps XMMWORD PTR [rax+rsi*1],xmm6
  95ce60:	mov    rax,QWORD PTR [rbx]
  95ce63:	mov    rsi,QWORD PTR [rsp+0x90]
  95ce6b:	vmovaps XMMWORD PTR [rax+rdi*1],xmm7
  95ce70:	mov    rax,QWORD PTR [rbx]
  95ce73:	vmovaps XMMWORD PTR [rax+rsi*1],xmm0
  95ce78:	add    QWORD PTR [rsp+0x68],0x4
  95ce7e:	mov    rdi,QWORD PTR [rsp+0x70]
  95ce83:	mov    rax,QWORD PTR [rsp+0x68]
  95ce88:	cmp    rax,QWORD PTR [rdi+0x40]
  95ce8c:	jb     95cb80 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x140>
  95ce92:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  95ce9e:	lea    rbp,[rax-0x360]
  95cea5:	mov    r12,rax
  95cea8:	mov    r13,QWORD PTR [rbp+0x188]
  95ceaf:	test   r13,r13
  95ceb2:	je     95ced2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x492>
  95ceb4:	movzx  eax,BYTE PTR [rbp+0x180]
  95cebb:	lea    ebx,[rax-0x1]
  95cebe:	mov    BYTE PTR [rbp+0x180],bl
  95cec4:	cmp    bl,0xf
  95cec7:	jbe    95cee4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4a4>
  95cec9:	add    WORD PTR [r13+0xfe2],0x1
  95ced2:	add    rsp,0xf8
  95ced9:	pop    rbx
  95ceda:	pop    rbp
  95cedb:	pop    r12
  95cedd:	pop    r13
  95cedf:	pop    r14
  95cee1:	pop    r15
  95cee3:	ret
  95cee4:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  95cee9:	movzx  r15d,bl
  95ceed:	lea    rdx,[r15+r15*2]
  95cef1:	mov    rcx,rax
  95cef4:	vmovq  xmm6,rax
  95cef9:	shl    rdx,0x3
  95cefd:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  95cf02:	vpinsrq xmm0,xmm6,rcx,0x1
  95cf08:	vpsubq xmm0,xmm0,XMMWORD PTR [r12+rdx*1-0x360]
  95cf12:	movzx  r12d,BYTE PTR [rbp+0x181]
  95cf1a:	cmp    r12b,bl
  95cf1d:	jae    95cf33 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4f3>
  95cf1f:	lea    eax,[r15-0x1]
  95cf23:	vmovq  rdx,xmm0
  95cf28:	cdqe
  95cf2a:	lea    rax,[rax+rax*2]
  95cf2e:	add    QWORD PTR [rbp+rax*8+0x8],rdx
  95cf33:	movzx  ebp,WORD PTR [r13+0xfe0]
  95cf3b:	cmp    bp,0x7f
  95cf3f:	je     95d09d <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x65d>
  95cf45:	vmovdqa XMMWORD PTR [rsp],xmm0
  95cf4a:	sub    ebx,r12d
  95cf4d:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  95cf59:	lea    rdx,[r15+r15*2]
  95cf5d:	vmovdqa xmm0,XMMWORD PTR [rsp]
  95cf62:	mov    DWORD PTR [rsp+0xe8],0x0
  95cf6d:	vmovdqu XMMWORD PTR [rsp+0xd8],xmm0
  95cf76:	lea    rdx,[rax+rdx*8-0x360]
  95cf7e:	mov    rax,QWORD PTR [rdx]
  95cf81:	mov    QWORD PTR [rsp+0xd0],rax
  95cf89:	movzx  eax,BYTE PTR [rdx+0x10]
  95cf8d:	lea    edx,[rbp+0x1]
  95cf90:	vmovdqa xmm7,XMMWORD PTR [rsp+0xd0]
  95cf99:	mov    WORD PTR [r13+0xfe0],dx
  95cfa1:	mov    ah,bl
  95cfa3:	mov    WORD PTR [rsp+0xec],ax
  95cfab:	movzx  eax,bp
  95cfae:	vmovdqu xmm6,XMMWORD PTR [rsp+0xde]
  95cfb7:	shl    rax,0x5
  95cfbb:	add    rax,r13
  95cfbe:	vmovdqu XMMWORD PTR [rax],xmm7
  95cfc2:	vmovdqu XMMWORD PTR [rax+0xe],xmm6
  95cfc7:	jmp    95ced2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x492>
  95cfcc:	shl    r8,0x4
  95cfd0:	add    rdi,r8
  95cfd3:	cmp    r9,rdi
  95cfd6:	je     95cb01 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc1>
  95cfdc:	mov    QWORD PTR [rbx+0x8],rdi
  95cfe0:	jmp    95cb01 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc1>
  95cfe5:	mov    eax,0x1
  95cfea:	lock xadd DWORD PTR [rip+0x2c9b6e],eax        # c26b60 <cpl::Profiling::registerRegion(char const*)::counter>
  95cff2:	mov    edx,0x1
  95cff7:	mov    r13d,0x1
  95cffd:	add    eax,0x2
  95d000:	cmp    eax,0xfe
  95d005:	ja     95d020 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x5e0>
  95d007:	lea    rdx,[rip+0x2c8012]        # c25020 <cpl::Profiling::regions>
  95d00e:	mov    ecx,eax
  95d010:	mov    r13d,eax
  95d013:	lea    rsi,[rip+0xa3298]        # a002b2 <_IO_stdin_used+0x182b2>
  95d01a:	mov    QWORD PTR [rdx+rcx*8],rsi
  95d01e:	mov    edx,eax
  95d020:	xor    eax,eax
  95d022:	lock cmpxchg BYTE PTR [rip+0x2d50be],dl        # c320e8 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  95d02a:	cmovne r13d,eax
  95d02e:	jmp    95ca79 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x39>
  95d033:	sub    rsi,rcx
  95d036:	mov    rdi,rbx
  95d039:	call   939b50 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  95d03e:	mov    rsi,QWORD PTR [rsp+0x70]
  95d043:	mov    rax,QWORD PTR [rsi+0x38]
  95d047:	mov    rdx,QWORD PTR [rsi+0x48]
  95d04b:	jmp    95cb01 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc1>
  95d050:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  95d055:	mov    BYTE PTR [rsp+0xe0],r13b
  95d05d:	mov    QWORD PTR [rsp+0xd0],rax
  95d065:	lea    rax,[r12+r12*2]
  95d069:	mov    QWORD PTR [rsp+0xd8],0x0
  95d075:	shl    rax,0x3
  95d079:	vmovdqa xmm7,XMMWORD PTR [rsp+0xd0]
  95d082:	vmovdqu XMMWORD PTR [rbp+rax*1-0x360],xmm7
  95d08b:	mov    BYTE PTR [rax+r14*1+0x10],r13b
  95d090:	movzx  r12d,BYTE PTR [r14+0x180]
  95d098:	jmp    95caab <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x6b>
  95d09d:	movzx  eax,WORD PTR [r13+0xfe2]
  95d0a5:	cmp    ax,0xffff
  95d0a9:	je     95ced2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x492>
  95d0af:	add    eax,0x1
  95d0b2:	mov    WORD PTR [r13+0xfe2],ax
  95d0ba:	jmp    95ced2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x492>
  95d0bf:	endbr64
  95d0c3:	mov    rbx,rax
  95d0c6:	vzeroupper
  95d0c9:	call   93c730 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  95d0ce:	mov    rdi,rbx
  95d0d1:	call   e1b10 <_Unwind_Resume@plt>

Disassembly of section .fini:
