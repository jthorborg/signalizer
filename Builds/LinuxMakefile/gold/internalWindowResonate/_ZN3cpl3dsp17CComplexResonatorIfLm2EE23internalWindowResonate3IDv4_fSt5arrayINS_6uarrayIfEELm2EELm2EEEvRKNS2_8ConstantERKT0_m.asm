; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

000000000095dc00 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long):
  95dc00:	endbr64
  95dc04:	push   r15
  95dc06:	mov    r15,rcx
  95dc09:	push   r14
  95dc0b:	push   r13
  95dc0d:	push   r12
  95dc0f:	push   rbp
  95dc10:	push   rbx
  95dc11:	mov    rbx,rdi
  95dc14:	sub    rsp,0xf8
  95dc1b:	mov    QWORD PTR [rsp+0x70],rsi
  95dc20:	mov    QWORD PTR [rsp+0xa0],rdx
  95dc28:	movzx  r13d,BYTE PTR [rip+0x2d54b8]        # c330e8 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  95dc30:	test   r13b,r13b
  95dc33:	je     95e1a5 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x5a5>
  95dc39:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  95dc45:	cmp    QWORD PTR [rax-0x1d8],0x0
  95dc4d:	mov    rbp,rax
  95dc50:	lea    r14,[rax-0x360]
  95dc57:	je     95dc82 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x82>
  95dc59:	movzx  r12d,BYTE PTR [r14+0x180]
  95dc61:	cmp    r12b,0xf
  95dc65:	jbe    95e210 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x610>
  95dc6b:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  95dc77:	add    r12d,0x1
  95dc7b:	mov    BYTE PTR [rax-0x1e0],r12b
  95dc82:	mov    rax,QWORD PTR [rsp+0x70]
  95dc87:	mov    r9,QWORD PTR [rbx+0x8]
  95dc8b:	mov    rdi,QWORD PTR [rbx]
  95dc8e:	mov    rdx,QWORD PTR [rax+0x48]
  95dc92:	mov    rax,QWORD PTR [rax+0x38]
  95dc96:	mov    rcx,r9
  95dc99:	sub    rcx,rdi
  95dc9c:	sar    rcx,0x2
  95dca0:	mov    r8,rax
  95dca3:	imul   r8,rdx
  95dca7:	lea    rsi,[r8*4+0x0]
  95dcaf:	cmp    rcx,rsi
  95dcb2:	jb     95e1f3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x5f3>
  95dcb8:	cmp    rsi,rcx
  95dcbb:	jb     95e18c <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x58c>
  95dcc1:	mov    rsi,QWORD PTR [rsp+0x70]
  95dcc6:	lea    rcx,[rdx+rdx*1]
  95dcca:	imul   rax,rcx
  95dcce:	cmp    QWORD PTR [rsi+0x40],0x0
  95dcd3:	je     95e052 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x452>
  95dcd9:	lea    r11,[rdx*4+0x0]
  95dce1:	xor    r10d,r10d
  95dce4:	mov    QWORD PTR [rsp+0x68],0x0
  95dced:	lea    rsi,[r11+rax*1]
  95dcf1:	shl    rax,0x2
  95dcf5:	mov    QWORD PTR [rsp+0xa8],r11
  95dcfd:	mov    QWORD PTR [rsp+0xc0],rax
  95dd05:	mov    rax,rsi
  95dd08:	lea    rdi,[rdx*8+0x0]
  95dd10:	sub    rax,rcx
  95dd13:	shl    rdx,0x4
  95dd17:	shl    rax,0x2
  95dd1b:	mov    QWORD PTR [rsp+0xb0],rdi
  95dd23:	mov    QWORD PTR [rsp+0xb8],rax
  95dd2b:	lea    rax,[rsi*4+0x0]
  95dd33:	mov    QWORD PTR [rsp+0xc8],rdx
  95dd3b:	mov    QWORD PTR [rsp+0x78],rax
  95dd40:	mov    rdi,QWORD PTR [rsp+0xb0]
  95dd48:	mov    rcx,QWORD PTR [rsp+0xa8]
  95dd50:	mov    rax,QWORD PTR [rsp+0x70]
  95dd55:	mov    rbp,QWORD PTR [rsp+0xc0]
  95dd5d:	add    rdi,r10
  95dd60:	mov    rdx,QWORD PTR [rsp+0xb8]
  95dd68:	lea    rsi,[rcx+rdi*1]
  95dd6c:	lea    r9,[rcx+rsi*1]
  95dd70:	mov    rcx,QWORD PTR [rsp+0xc8]
  95dd78:	mov    rax,QWORD PTR [rax]
  95dd7b:	lea    r13,[rbp+r10*1+0x0]
  95dd80:	lea    r12,[rbp+r11*1+0x0]
  95dd85:	mov    rbp,rdx
  95dd88:	add    rdx,r11
  95dd8b:	add    rbp,r10
  95dd8e:	mov    QWORD PTR [rsp+0x80],rdx
  95dd96:	lea    r8,[rcx+r11*1]
  95dd9a:	mov    rcx,QWORD PTR [rsp+0xa0]
  95dda2:	vmovaps xmm6,XMMWORD PTR [rax+r10*1]
  95dda8:	vmovaps xmm7,XMMWORD PTR [rax+r11*1]
  95ddae:	vmovaps xmm9,XMMWORD PTR [rax+rdi*1]
  95ddb3:	vmovaps xmm10,XMMWORD PTR [rax+rsi*1]
  95ddb8:	vmovaps xmm11,XMMWORD PTR [rax+r9*1]
  95ddbe:	vmovaps xmm8,XMMWORD PTR [rax+r8*1]
  95ddc4:	mov    rax,QWORD PTR [rcx]
  95ddc7:	mov    rcx,QWORD PTR [rcx+0x10]
  95ddcb:	mov    QWORD PTR [rsp+0x98],rax
  95ddd3:	mov    rax,QWORD PTR [rbx]
  95ddd6:	vmovaps xmm0,XMMWORD PTR [rax+r8*1]
  95dddc:	lea    r14,[rax+r10*1]
  95dde0:	vmovaps xmm5,XMMWORD PTR [r14]
  95dde5:	vmovaps XMMWORD PTR [rsp+0x30],xmm0
  95ddeb:	vmovaps xmm0,XMMWORD PTR [rax+r13*1]
  95ddf1:	vmovaps xmm4,XMMWORD PTR [rax+r11*1]
  95ddf7:	vmovaps XMMWORD PTR [rsp+0x10],xmm0
  95ddfd:	vmovaps xmm0,XMMWORD PTR [rax+r12*1]
  95de03:	vmovaps xmm3,XMMWORD PTR [rax+rdi*1]
  95de08:	vmovaps XMMWORD PTR [rsp+0x50],xmm0
  95de0e:	vmovaps xmm0,XMMWORD PTR [rax+rbp*1]
  95de13:	vmovaps xmm2,XMMWORD PTR [rax+rsi*1]
  95de18:	vmovaps xmm1,XMMWORD PTR [rax+r9*1]
  95de1e:	vmovaps XMMWORD PTR [rsp],xmm0
  95de23:	vmovaps xmm0,XMMWORD PTR [rax+rdx*1]
  95de28:	mov    rdx,QWORD PTR [rsp+0x78]
  95de2d:	vmovaps XMMWORD PTR [rsp+0x40],xmm0
  95de33:	add    rdx,r10
  95de36:	vmovaps xmm0,XMMWORD PTR [rax+rdx*1]
  95de3b:	mov    QWORD PTR [rsp+0x88],rdx
  95de43:	mov    rdx,QWORD PTR [rsp+0x78]
  95de48:	vmovaps XMMWORD PTR [rsp+0x20],xmm0
  95de4e:	add    rdx,r11
  95de51:	vmovaps xmm0,XMMWORD PTR [rax+rdx*1]
  95de56:	mov    QWORD PTR [rsp+0x90],rdx
  95de5e:	test   r15,r15
  95de61:	je     95df93 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x393>
  95de67:	mov    rdx,QWORD PTR [rsp+0x98]
  95de6f:	xor    eax,eax
  95de71:	nop    DWORD PTR [rax+0x0]
  95de78:	vmulps xmm12,xmm7,xmm4
  95de7c:	vbroadcastss xmm13,DWORD PTR [rdx+rax*4]
  95de82:	vmulps xmm14,xmm5,xmm7
  95de86:	vmovaps xmm15,XMMWORD PTR [rsp+0x30]
  95de8c:	vmulps xmm5,xmm5,xmm6
  95de90:	vmulps xmm4,xmm6,xmm4
  95de94:	vsubps xmm12,xmm13,xmm12
  95de99:	vaddps xmm5,xmm5,xmm12
  95de9e:	vmulps xmm12,xmm2,xmm10
  95dea3:	vaddps xmm4,xmm14,xmm4
  95dea7:	vmulps xmm14,xmm10,xmm3
  95deab:	vmulps xmm3,xmm9,xmm3
  95deaf:	vmulps xmm2,xmm2,xmm9
  95deb4:	vsubps xmm12,xmm13,xmm12
  95deb9:	vaddps xmm3,xmm3,xmm12
  95debe:	vmulps xmm12,xmm8,xmm15
  95dec3:	vaddps xmm2,xmm14,xmm2
  95dec7:	vmulps xmm14,xmm8,xmm1
  95decb:	vmulps xmm1,xmm11,xmm1
  95decf:	vsubps xmm12,xmm13,xmm12
  95ded4:	vmulps xmm13,xmm11,xmm15
  95ded9:	vmulps xmm15,xmm7,XMMWORD PTR [rsp+0x10]
  95dedf:	vaddps xmm1,xmm1,xmm12
  95dee4:	vbroadcastss xmm12,DWORD PTR [rcx+rax*4]
  95deea:	add    rax,0x1
  95deee:	vaddps xmm13,xmm14,xmm13
  95def3:	vmovaps xmm14,XMMWORD PTR [rsp+0x50]
  95def9:	vmovaps XMMWORD PTR [rsp+0x30],xmm13
  95deff:	vmulps xmm13,xmm7,xmm14
  95df04:	vmulps xmm14,xmm6,xmm14
  95df09:	vsubps xmm13,xmm12,xmm13
  95df0e:	vaddps xmm14,xmm15,xmm14
  95df13:	vmulps xmm15,xmm10,XMMWORD PTR [rsp]
  95df18:	vmovaps XMMWORD PTR [rsp+0x50],xmm14
  95df1e:	vmulps xmm14,xmm6,XMMWORD PTR [rsp+0x10]
  95df24:	vaddps xmm13,xmm14,xmm13
  95df29:	vmovaps xmm14,XMMWORD PTR [rsp+0x40]
  95df2f:	vmovaps XMMWORD PTR [rsp+0x10],xmm13
  95df35:	vmulps xmm13,xmm10,xmm14
  95df3a:	vmulps xmm14,xmm9,xmm14
  95df3f:	vsubps xmm13,xmm12,xmm13
  95df44:	vaddps xmm14,xmm15,xmm14
  95df49:	vmovaps xmm15,XMMWORD PTR [rsp+0x20]
  95df4f:	vmovaps XMMWORD PTR [rsp+0x40],xmm14
  95df55:	vmulps xmm14,xmm9,XMMWORD PTR [rsp]
  95df5a:	vaddps xmm13,xmm14,xmm13
  95df5f:	vmovaps XMMWORD PTR [rsp],xmm13
  95df64:	vmulps xmm13,xmm8,xmm0
  95df68:	vmulps xmm0,xmm11,xmm0
  95df6c:	vsubps xmm12,xmm12,xmm13
  95df71:	vmulps xmm13,xmm8,xmm15
  95df76:	vaddps xmm0,xmm13,xmm0
  95df7a:	vmulps xmm13,xmm11,xmm15
  95df7f:	vaddps xmm12,xmm13,xmm12
  95df84:	vmovaps XMMWORD PTR [rsp+0x20],xmm12
  95df8a:	cmp    r15,rax
  95df8d:	jne    95de78 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x278>
  95df93:	vmovaps XMMWORD PTR [r14],xmm5
  95df98:	mov    rax,QWORD PTR [rbx]
  95df9b:	add    r10,0x10
  95df9f:	vmovaps xmm6,XMMWORD PTR [rsp+0x30]
  95dfa5:	vmovaps xmm7,XMMWORD PTR [rsp+0x10]
  95dfab:	vmovaps XMMWORD PTR [rax+r11*1],xmm4
  95dfb1:	mov    rax,QWORD PTR [rbx]
  95dfb4:	add    r11,0x10
  95dfb8:	vmovaps XMMWORD PTR [rax+rdi*1],xmm3
  95dfbd:	mov    rax,QWORD PTR [rbx]
  95dfc0:	mov    rdi,QWORD PTR [rsp+0x88]
  95dfc8:	vmovaps XMMWORD PTR [rax+rsi*1],xmm2
  95dfcd:	mov    rax,QWORD PTR [rbx]
  95dfd0:	mov    rsi,QWORD PTR [rsp+0x80]
  95dfd8:	vmovaps XMMWORD PTR [rax+r9*1],xmm1
  95dfde:	mov    rax,QWORD PTR [rbx]
  95dfe1:	vmovaps XMMWORD PTR [rax+r8*1],xmm6
  95dfe7:	mov    rax,QWORD PTR [rbx]
  95dfea:	vmovaps xmm6,XMMWORD PTR [rsp+0x50]
  95dff0:	vmovaps XMMWORD PTR [rax+r13*1],xmm7
  95dff6:	mov    rax,QWORD PTR [rbx]
  95dff9:	vmovaps xmm7,XMMWORD PTR [rsp]
  95dffe:	vmovaps XMMWORD PTR [rax+r12*1],xmm6
  95e004:	mov    rax,QWORD PTR [rbx]
  95e007:	vmovaps xmm6,XMMWORD PTR [rsp+0x40]
  95e00d:	vmovaps XMMWORD PTR [rax+rbp*1],xmm7
  95e012:	mov    rax,QWORD PTR [rbx]
  95e015:	vmovaps xmm7,XMMWORD PTR [rsp+0x20]
  95e01b:	vmovaps XMMWORD PTR [rax+rsi*1],xmm6
  95e020:	mov    rax,QWORD PTR [rbx]
  95e023:	mov    rsi,QWORD PTR [rsp+0x90]
  95e02b:	vmovaps XMMWORD PTR [rax+rdi*1],xmm7
  95e030:	mov    rax,QWORD PTR [rbx]
  95e033:	vmovaps XMMWORD PTR [rax+rsi*1],xmm0
  95e038:	add    QWORD PTR [rsp+0x68],0x4
  95e03e:	mov    rdi,QWORD PTR [rsp+0x70]
  95e043:	mov    rax,QWORD PTR [rsp+0x68]
  95e048:	cmp    rax,QWORD PTR [rdi+0x40]
  95e04c:	jb     95dd40 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x140>
  95e052:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  95e05e:	lea    rbp,[rax-0x360]
  95e065:	mov    r12,rax
  95e068:	mov    r13,QWORD PTR [rbp+0x188]
  95e06f:	test   r13,r13
  95e072:	je     95e092 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x492>
  95e074:	movzx  eax,BYTE PTR [rbp+0x180]
  95e07b:	lea    ebx,[rax-0x1]
  95e07e:	mov    BYTE PTR [rbp+0x180],bl
  95e084:	cmp    bl,0xf
  95e087:	jbe    95e0a4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4a4>
  95e089:	add    WORD PTR [r13+0xfe2],0x1
  95e092:	add    rsp,0xf8
  95e099:	pop    rbx
  95e09a:	pop    rbp
  95e09b:	pop    r12
  95e09d:	pop    r13
  95e09f:	pop    r14
  95e0a1:	pop    r15
  95e0a3:	ret
  95e0a4:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  95e0a9:	movzx  r15d,bl
  95e0ad:	lea    rdx,[r15+r15*2]
  95e0b1:	mov    rcx,rax
  95e0b4:	vmovq  xmm6,rax
  95e0b9:	shl    rdx,0x3
  95e0bd:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  95e0c2:	vpinsrq xmm0,xmm6,rcx,0x1
  95e0c8:	vpsubq xmm0,xmm0,XMMWORD PTR [r12+rdx*1-0x360]
  95e0d2:	movzx  r12d,BYTE PTR [rbp+0x181]
  95e0da:	cmp    r12b,bl
  95e0dd:	jae    95e0f3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x4f3>
  95e0df:	lea    eax,[r15-0x1]
  95e0e3:	vmovq  rdx,xmm0
  95e0e8:	cdqe
  95e0ea:	lea    rax,[rax+rax*2]
  95e0ee:	add    QWORD PTR [rbp+rax*8+0x8],rdx
  95e0f3:	movzx  ebp,WORD PTR [r13+0xfe0]
  95e0fb:	cmp    bp,0x7f
  95e0ff:	je     95e25d <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x65d>
  95e105:	vmovdqa XMMWORD PTR [rsp],xmm0
  95e10a:	sub    ebx,r12d
  95e10d:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  95e119:	lea    rdx,[r15+r15*2]
  95e11d:	vmovdqa xmm0,XMMWORD PTR [rsp]
  95e122:	mov    DWORD PTR [rsp+0xe8],0x0
  95e12d:	vmovdqu XMMWORD PTR [rsp+0xd8],xmm0
  95e136:	lea    rdx,[rax+rdx*8-0x360]
  95e13e:	mov    rax,QWORD PTR [rdx]
  95e141:	mov    QWORD PTR [rsp+0xd0],rax
  95e149:	movzx  eax,BYTE PTR [rdx+0x10]
  95e14d:	lea    edx,[rbp+0x1]
  95e150:	vmovdqa xmm7,XMMWORD PTR [rsp+0xd0]
  95e159:	mov    WORD PTR [r13+0xfe0],dx
  95e161:	mov    ah,bl
  95e163:	mov    WORD PTR [rsp+0xec],ax
  95e16b:	movzx  eax,bp
  95e16e:	vmovdqu xmm6,XMMWORD PTR [rsp+0xde]
  95e177:	shl    rax,0x5
  95e17b:	add    rax,r13
  95e17e:	vmovdqu XMMWORD PTR [rax],xmm7
  95e182:	vmovdqu XMMWORD PTR [rax+0xe],xmm6
  95e187:	jmp    95e092 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x492>
  95e18c:	shl    r8,0x4
  95e190:	add    rdi,r8
  95e193:	cmp    r9,rdi
  95e196:	je     95dcc1 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc1>
  95e19c:	mov    QWORD PTR [rbx+0x8],rdi
  95e1a0:	jmp    95dcc1 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc1>
  95e1a5:	mov    eax,0x1
  95e1aa:	lock xadd DWORD PTR [rip+0x2c99ae],eax        # c27b60 <cpl::Profiling::registerRegion(char const*)::counter>
  95e1b2:	mov    edx,0x1
  95e1b7:	mov    r13d,0x1
  95e1bd:	add    eax,0x2
  95e1c0:	cmp    eax,0xfe
  95e1c5:	ja     95e1e0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x5e0>
  95e1c7:	lea    rdx,[rip+0x2c7e52]        # c26020 <cpl::Profiling::regions>
  95e1ce:	mov    ecx,eax
  95e1d0:	mov    r13d,eax
  95e1d3:	lea    rsi,[rip+0xa30d8]        # a012b2 <_IO_stdin_used+0x182b2>
  95e1da:	mov    QWORD PTR [rdx+rcx*8],rsi
  95e1de:	mov    edx,eax
  95e1e0:	xor    eax,eax
  95e1e2:	lock cmpxchg BYTE PTR [rip+0x2d4efe],dl        # c330e8 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)::profilerCached49>
  95e1ea:	cmovne r13d,eax
  95e1ee:	jmp    95dc39 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x39>
  95e1f3:	sub    rsi,rcx
  95e1f6:	mov    rdi,rbx
  95e1f9:	call   939bb0 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  95e1fe:	mov    rsi,QWORD PTR [rsp+0x70]
  95e203:	mov    rax,QWORD PTR [rsi+0x38]
  95e207:	mov    rdx,QWORD PTR [rsi+0x48]
  95e20b:	jmp    95dcc1 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0xc1>
  95e210:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  95e215:	mov    BYTE PTR [rsp+0xe0],r13b
  95e21d:	mov    QWORD PTR [rsp+0xd0],rax
  95e225:	lea    rax,[r12+r12*2]
  95e229:	mov    QWORD PTR [rsp+0xd8],0x0
  95e235:	shl    rax,0x3
  95e239:	vmovdqa xmm7,XMMWORD PTR [rsp+0xd0]
  95e242:	vmovdqu XMMWORD PTR [rbp+rax*1-0x360],xmm7
  95e24b:	mov    BYTE PTR [rax+r14*1+0x10],r13b
  95e250:	movzx  r12d,BYTE PTR [r14+0x180]
  95e258:	jmp    95dc6b <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x6b>
  95e25d:	movzx  eax,WORD PTR [r13+0xfe2]
  95e265:	cmp    ax,0xffff
  95e269:	je     95e092 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x492>
  95e26f:	add    eax,0x1
  95e272:	mov    WORD PTR [r13+0xfe2],ax
  95e27a:	jmp    95e092 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 2ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 2ul> const&, unsigned long)+0x492>
  95e27f:	endbr64
  95e283:	mov    rbx,rax
  95e286:	vzeroupper
  95e289:	call   93c790 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  95e28e:	mov    rdi,rbx
  95e291:	call   e1b10 <_Unwind_Resume@plt>

Disassembly of section .fini:
