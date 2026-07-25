; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000959a20 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long):
  959a20:	endbr64
  959a24:	push   r15
  959a26:	mov    r15,rcx
  959a29:	push   r14
  959a2b:	push   r13
  959a2d:	push   r12
  959a2f:	push   rbp
  959a30:	push   rbx
  959a31:	mov    rbx,rdi
  959a34:	sub    rsp,0xf8
  959a3b:	mov    QWORD PTR [rsp+0x70],rsi
  959a40:	mov    QWORD PTR [rsp+0xa0],rdx
  959a48:	movzx  r13d,BYTE PTR [rip+0x2d96a2]        # c330f2 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  959a50:	test   r13b,r13b
  959a53:	je     959fc5 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x5a5>
  959a59:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  959a65:	cmp    QWORD PTR [rax-0x1d8],0x0
  959a6d:	mov    rbp,rax
  959a70:	lea    r14,[rax-0x360]
  959a77:	je     959aa2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x82>
  959a79:	movzx  r12d,BYTE PTR [r14+0x180]
  959a81:	cmp    r12b,0xf
  959a85:	jbe    95a030 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x610>
  959a8b:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  959a97:	add    r12d,0x1
  959a9b:	mov    BYTE PTR [rax-0x1e0],r12b
  959aa2:	mov    rax,QWORD PTR [rsp+0x70]
  959aa7:	mov    r9,QWORD PTR [rbx+0x8]
  959aab:	mov    rdi,QWORD PTR [rbx]
  959aae:	mov    rdx,QWORD PTR [rax+0x48]
  959ab2:	mov    rax,QWORD PTR [rax+0x38]
  959ab6:	mov    rcx,r9
  959ab9:	sub    rcx,rdi
  959abc:	sar    rcx,0x2
  959ac0:	mov    r8,rax
  959ac3:	imul   r8,rdx
  959ac7:	lea    rsi,[r8*4+0x0]
  959acf:	cmp    rcx,rsi
  959ad2:	jb     95a013 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x5f3>
  959ad8:	cmp    rsi,rcx
  959adb:	jb     959fac <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x58c>
  959ae1:	mov    rsi,QWORD PTR [rsp+0x70]
  959ae6:	lea    rcx,[rdx+rdx*1]
  959aea:	imul   rax,rcx
  959aee:	cmp    QWORD PTR [rsi+0x40],0x0
  959af3:	je     959e72 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x452>
  959af9:	lea    r11,[rdx*4+0x0]
  959b01:	xor    r10d,r10d
  959b04:	mov    QWORD PTR [rsp+0x68],0x0
  959b0d:	lea    rsi,[r11+rax*1]
  959b11:	shl    rax,0x2
  959b15:	mov    QWORD PTR [rsp+0xa8],r11
  959b1d:	mov    QWORD PTR [rsp+0xc0],rax
  959b25:	mov    rax,rsi
  959b28:	lea    rdi,[rdx*8+0x0]
  959b30:	sub    rax,rcx
  959b33:	shl    rdx,0x4
  959b37:	shl    rax,0x2
  959b3b:	mov    QWORD PTR [rsp+0xb0],rdi
  959b43:	mov    QWORD PTR [rsp+0xb8],rax
  959b4b:	lea    rax,[rsi*4+0x0]
  959b53:	mov    QWORD PTR [rsp+0xc8],rdx
  959b5b:	mov    QWORD PTR [rsp+0x78],rax
  959b60:	mov    rdi,QWORD PTR [rsp+0xb0]
  959b68:	mov    rcx,QWORD PTR [rsp+0xa8]
  959b70:	mov    rax,QWORD PTR [rsp+0x70]
  959b75:	mov    rbp,QWORD PTR [rsp+0xc0]
  959b7d:	add    rdi,r10
  959b80:	mov    rdx,QWORD PTR [rsp+0xb8]
  959b88:	lea    rsi,[rcx+rdi*1]
  959b8c:	lea    r9,[rcx+rsi*1]
  959b90:	mov    rcx,QWORD PTR [rsp+0xc8]
  959b98:	mov    rax,QWORD PTR [rax]
  959b9b:	lea    r13,[rbp+r10*1+0x0]
  959ba0:	lea    r12,[rbp+r11*1+0x0]
  959ba5:	mov    rbp,rdx
  959ba8:	add    rdx,r11
  959bab:	add    rbp,r10
  959bae:	mov    QWORD PTR [rsp+0x80],rdx
  959bb6:	lea    r8,[rcx+r11*1]
  959bba:	mov    rcx,QWORD PTR [rsp+0xa0]
  959bc2:	vmovaps xmm6,XMMWORD PTR [rax+r10*1]
  959bc8:	vmovaps xmm7,XMMWORD PTR [rax+r11*1]
  959bce:	vmovaps xmm9,XMMWORD PTR [rax+rdi*1]
  959bd3:	vmovaps xmm10,XMMWORD PTR [rax+rsi*1]
  959bd8:	vmovaps xmm11,XMMWORD PTR [rax+r9*1]
  959bde:	vmovaps xmm8,XMMWORD PTR [rax+r8*1]
  959be4:	mov    rax,QWORD PTR [rcx]
  959be7:	mov    rcx,QWORD PTR [rcx+0x10]
  959beb:	mov    QWORD PTR [rsp+0x98],rax
  959bf3:	mov    rax,QWORD PTR [rbx]
  959bf6:	vmovaps xmm0,XMMWORD PTR [rax+r8*1]
  959bfc:	lea    r14,[rax+r10*1]
  959c00:	vmovaps xmm5,XMMWORD PTR [r14]
  959c05:	vmovaps XMMWORD PTR [rsp+0x30],xmm0
  959c0b:	vmovaps xmm0,XMMWORD PTR [rax+r13*1]
  959c11:	vmovaps xmm4,XMMWORD PTR [rax+r11*1]
  959c17:	vmovaps XMMWORD PTR [rsp+0x10],xmm0
  959c1d:	vmovaps xmm0,XMMWORD PTR [rax+r12*1]
  959c23:	vmovaps xmm3,XMMWORD PTR [rax+rdi*1]
  959c28:	vmovaps XMMWORD PTR [rsp+0x50],xmm0
  959c2e:	vmovaps xmm0,XMMWORD PTR [rax+rbp*1]
  959c33:	vmovaps xmm2,XMMWORD PTR [rax+rsi*1]
  959c38:	vmovaps xmm1,XMMWORD PTR [rax+r9*1]
  959c3e:	vmovaps XMMWORD PTR [rsp],xmm0
  959c43:	vmovaps xmm0,XMMWORD PTR [rax+rdx*1]
  959c48:	mov    rdx,QWORD PTR [rsp+0x78]
  959c4d:	vmovaps XMMWORD PTR [rsp+0x40],xmm0
  959c53:	add    rdx,r10
  959c56:	vmovaps xmm0,XMMWORD PTR [rax+rdx*1]
  959c5b:	mov    QWORD PTR [rsp+0x88],rdx
  959c63:	mov    rdx,QWORD PTR [rsp+0x78]
  959c68:	vmovaps XMMWORD PTR [rsp+0x20],xmm0
  959c6e:	add    rdx,r11
  959c71:	vmovaps xmm0,XMMWORD PTR [rax+rdx*1]
  959c76:	mov    QWORD PTR [rsp+0x90],rdx
  959c7e:	test   r15,r15
  959c81:	je     959db3 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x393>
  959c87:	mov    rdx,QWORD PTR [rsp+0x98]
  959c8f:	xor    eax,eax
  959c91:	nop    DWORD PTR [rax+0x0]
  959c98:	vmulps xmm12,xmm7,xmm4
  959c9c:	vbroadcastss xmm13,DWORD PTR [rdx+rax*4]
  959ca2:	vmulps xmm14,xmm5,xmm7
  959ca6:	vmovaps xmm15,XMMWORD PTR [rsp+0x30]
  959cac:	vmulps xmm5,xmm5,xmm6
  959cb0:	vmulps xmm4,xmm6,xmm4
  959cb4:	vsubps xmm12,xmm13,xmm12
  959cb9:	vaddps xmm5,xmm5,xmm12
  959cbe:	vmulps xmm12,xmm2,xmm10
  959cc3:	vaddps xmm4,xmm14,xmm4
  959cc7:	vmulps xmm14,xmm10,xmm3
  959ccb:	vmulps xmm3,xmm9,xmm3
  959ccf:	vmulps xmm2,xmm2,xmm9
  959cd4:	vsubps xmm12,xmm13,xmm12
  959cd9:	vaddps xmm3,xmm3,xmm12
  959cde:	vmulps xmm12,xmm8,xmm15
  959ce3:	vaddps xmm2,xmm14,xmm2
  959ce7:	vmulps xmm14,xmm8,xmm1
  959ceb:	vmulps xmm1,xmm11,xmm1
  959cef:	vsubps xmm12,xmm13,xmm12
  959cf4:	vmulps xmm13,xmm11,xmm15
  959cf9:	vmulps xmm15,xmm7,XMMWORD PTR [rsp+0x10]
  959cff:	vaddps xmm1,xmm1,xmm12
  959d04:	vbroadcastss xmm12,DWORD PTR [rcx+rax*4]
  959d0a:	add    rax,0x1
  959d0e:	vaddps xmm13,xmm14,xmm13
  959d13:	vmovaps xmm14,XMMWORD PTR [rsp+0x50]
  959d19:	vmovaps XMMWORD PTR [rsp+0x30],xmm13
  959d1f:	vmulps xmm13,xmm7,xmm14
  959d24:	vmulps xmm14,xmm6,xmm14
  959d29:	vsubps xmm13,xmm12,xmm13
  959d2e:	vaddps xmm14,xmm15,xmm14
  959d33:	vmulps xmm15,xmm10,XMMWORD PTR [rsp]
  959d38:	vmovaps XMMWORD PTR [rsp+0x50],xmm14
  959d3e:	vmulps xmm14,xmm6,XMMWORD PTR [rsp+0x10]
  959d44:	vaddps xmm13,xmm14,xmm13
  959d49:	vmovaps xmm14,XMMWORD PTR [rsp+0x40]
  959d4f:	vmovaps XMMWORD PTR [rsp+0x10],xmm13
  959d55:	vmulps xmm13,xmm10,xmm14
  959d5a:	vmulps xmm14,xmm9,xmm14
  959d5f:	vsubps xmm13,xmm12,xmm13
  959d64:	vaddps xmm14,xmm15,xmm14
  959d69:	vmovaps xmm15,XMMWORD PTR [rsp+0x20]
  959d6f:	vmovaps XMMWORD PTR [rsp+0x40],xmm14
  959d75:	vmulps xmm14,xmm9,XMMWORD PTR [rsp]
  959d7a:	vaddps xmm13,xmm14,xmm13
  959d7f:	vmovaps XMMWORD PTR [rsp],xmm13
  959d84:	vmulps xmm13,xmm8,xmm0
  959d88:	vmulps xmm0,xmm11,xmm0
  959d8c:	vsubps xmm12,xmm12,xmm13
  959d91:	vmulps xmm13,xmm8,xmm15
  959d96:	vaddps xmm0,xmm13,xmm0
  959d9a:	vmulps xmm13,xmm11,xmm15
  959d9f:	vaddps xmm12,xmm13,xmm12
  959da4:	vmovaps XMMWORD PTR [rsp+0x20],xmm12
  959daa:	cmp    r15,rax
  959dad:	jne    959c98 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x278>
  959db3:	vmovaps XMMWORD PTR [r14],xmm5
  959db8:	mov    rax,QWORD PTR [rbx]
  959dbb:	add    r10,0x10
  959dbf:	vmovaps xmm6,XMMWORD PTR [rsp+0x30]
  959dc5:	vmovaps xmm7,XMMWORD PTR [rsp+0x10]
  959dcb:	vmovaps XMMWORD PTR [rax+r11*1],xmm4
  959dd1:	mov    rax,QWORD PTR [rbx]
  959dd4:	add    r11,0x10
  959dd8:	vmovaps XMMWORD PTR [rax+rdi*1],xmm3
  959ddd:	mov    rax,QWORD PTR [rbx]
  959de0:	mov    rdi,QWORD PTR [rsp+0x88]
  959de8:	vmovaps XMMWORD PTR [rax+rsi*1],xmm2
  959ded:	mov    rax,QWORD PTR [rbx]
  959df0:	mov    rsi,QWORD PTR [rsp+0x80]
  959df8:	vmovaps XMMWORD PTR [rax+r9*1],xmm1
  959dfe:	mov    rax,QWORD PTR [rbx]
  959e01:	vmovaps XMMWORD PTR [rax+r8*1],xmm6
  959e07:	mov    rax,QWORD PTR [rbx]
  959e0a:	vmovaps xmm6,XMMWORD PTR [rsp+0x50]
  959e10:	vmovaps XMMWORD PTR [rax+r13*1],xmm7
  959e16:	mov    rax,QWORD PTR [rbx]
  959e19:	vmovaps xmm7,XMMWORD PTR [rsp]
  959e1e:	vmovaps XMMWORD PTR [rax+r12*1],xmm6
  959e24:	mov    rax,QWORD PTR [rbx]
  959e27:	vmovaps xmm6,XMMWORD PTR [rsp+0x40]
  959e2d:	vmovaps XMMWORD PTR [rax+rbp*1],xmm7
  959e32:	mov    rax,QWORD PTR [rbx]
  959e35:	vmovaps xmm7,XMMWORD PTR [rsp+0x20]
  959e3b:	vmovaps XMMWORD PTR [rax+rsi*1],xmm6
  959e40:	mov    rax,QWORD PTR [rbx]
  959e43:	mov    rsi,QWORD PTR [rsp+0x90]
  959e4b:	vmovaps XMMWORD PTR [rax+rdi*1],xmm7
  959e50:	mov    rax,QWORD PTR [rbx]
  959e53:	vmovaps XMMWORD PTR [rax+rsi*1],xmm0
  959e58:	add    QWORD PTR [rsp+0x68],0x4
  959e5e:	mov    rdi,QWORD PTR [rsp+0x70]
  959e63:	mov    rax,QWORD PTR [rsp+0x68]
  959e68:	cmp    rax,QWORD PTR [rdi+0x40]
  959e6c:	jb     959b60 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x140>
  959e72:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  959e7e:	lea    rbp,[rax-0x360]
  959e85:	mov    r12,rax
  959e88:	mov    r13,QWORD PTR [rbp+0x188]
  959e8f:	test   r13,r13
  959e92:	je     959eb2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x492>
  959e94:	movzx  eax,BYTE PTR [rbp+0x180]
  959e9b:	lea    ebx,[rax-0x1]
  959e9e:	mov    BYTE PTR [rbp+0x180],bl
  959ea4:	cmp    bl,0xf
  959ea7:	jbe    959ec4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4a4>
  959ea9:	add    WORD PTR [r13+0xfe2],0x1
  959eb2:	add    rsp,0xf8
  959eb9:	pop    rbx
  959eba:	pop    rbp
  959ebb:	pop    r12
  959ebd:	pop    r13
  959ebf:	pop    r14
  959ec1:	pop    r15
  959ec3:	ret
  959ec4:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  959ec9:	movzx  r15d,bl
  959ecd:	lea    rdx,[r15+r15*2]
  959ed1:	mov    rcx,rax
  959ed4:	vmovq  xmm6,rax
  959ed9:	shl    rdx,0x3
  959edd:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  959ee2:	vpinsrq xmm0,xmm6,rcx,0x1
  959ee8:	vpsubq xmm0,xmm0,XMMWORD PTR [r12+rdx*1-0x360]
  959ef2:	movzx  r12d,BYTE PTR [rbp+0x181]
  959efa:	cmp    r12b,bl
  959efd:	jae    959f13 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4f3>
  959eff:	lea    eax,[r15-0x1]
  959f03:	vmovq  rdx,xmm0
  959f08:	cdqe
  959f0a:	lea    rax,[rax+rax*2]
  959f0e:	add    QWORD PTR [rbp+rax*8+0x8],rdx
  959f13:	movzx  ebp,WORD PTR [r13+0xfe0]
  959f1b:	cmp    bp,0x7f
  959f1f:	je     95a07d <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x65d>
  959f25:	vmovdqa XMMWORD PTR [rsp],xmm0
  959f2a:	sub    ebx,r12d
  959f2d:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  959f39:	lea    rdx,[r15+r15*2]
  959f3d:	vmovdqa xmm0,XMMWORD PTR [rsp]
  959f42:	mov    DWORD PTR [rsp+0xe8],0x0
  959f4d:	vmovdqu XMMWORD PTR [rsp+0xd8],xmm0
  959f56:	lea    rdx,[rax+rdx*8-0x360]
  959f5e:	mov    rax,QWORD PTR [rdx]
  959f61:	mov    QWORD PTR [rsp+0xd0],rax
  959f69:	movzx  eax,BYTE PTR [rdx+0x10]
  959f6d:	lea    edx,[rbp+0x1]
  959f70:	vmovdqa xmm7,XMMWORD PTR [rsp+0xd0]
  959f79:	mov    WORD PTR [r13+0xfe0],dx
  959f81:	mov    ah,bl
  959f83:	mov    WORD PTR [rsp+0xec],ax
  959f8b:	movzx  eax,bp
  959f8e:	vmovdqu xmm6,XMMWORD PTR [rsp+0xde]
  959f97:	shl    rax,0x5
  959f9b:	add    rax,r13
  959f9e:	vmovdqu XMMWORD PTR [rax],xmm7
  959fa2:	vmovdqu XMMWORD PTR [rax+0xe],xmm6
  959fa7:	jmp    959eb2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x492>
  959fac:	shl    r8,0x4
  959fb0:	add    rdi,r8
  959fb3:	cmp    r9,rdi
  959fb6:	je     959ae1 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc1>
  959fbc:	mov    QWORD PTR [rbx+0x8],rdi
  959fc0:	jmp    959ae1 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc1>
  959fc5:	mov    eax,0x1
  959fca:	lock xadd DWORD PTR [rip+0x2cdb8e],eax        # c27b60 <cpl::Profiling::registerRegion(char const*)::counter>
  959fd2:	mov    edx,0x1
  959fd7:	mov    r13d,0x1
  959fdd:	add    eax,0x2
  959fe0:	cmp    eax,0xfe
  959fe5:	ja     95a000 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x5e0>
  959fe7:	lea    rdx,[rip+0x2cc032]        # c26020 <cpl::Profiling::regions>
  959fee:	mov    ecx,eax
  959ff0:	mov    r13d,eax
  959ff3:	lea    rsi,[rip+0xa72b8]        # a012b2 <_IO_stdin_used+0x182b2>
  959ffa:	mov    QWORD PTR [rdx+rcx*8],rsi
  959ffe:	mov    edx,eax
  95a000:	xor    eax,eax
  95a002:	lock cmpxchg BYTE PTR [rip+0x2d90e8],dl        # c330f2 <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  95a00a:	cmovne r13d,eax
  95a00e:	jmp    959a59 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x39>
  95a013:	sub    rsi,rcx
  95a016:	mov    rdi,rbx
  95a019:	call   939bb0 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  95a01e:	mov    rsi,QWORD PTR [rsp+0x70]
  95a023:	mov    rax,QWORD PTR [rsi+0x38]
  95a027:	mov    rdx,QWORD PTR [rsi+0x48]
  95a02b:	jmp    959ae1 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xc1>
  95a030:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  95a035:	mov    BYTE PTR [rsp+0xe0],r13b
  95a03d:	mov    QWORD PTR [rsp+0xd0],rax
  95a045:	lea    rax,[r12+r12*2]
  95a049:	mov    QWORD PTR [rsp+0xd8],0x0
  95a055:	shl    rax,0x3
  95a059:	vmovdqa xmm7,XMMWORD PTR [rsp+0xd0]
  95a062:	vmovdqu XMMWORD PTR [rbp+rax*1-0x360],xmm7
  95a06b:	mov    BYTE PTR [rax+r14*1+0x10],r13b
  95a070:	movzx  r12d,BYTE PTR [r14+0x180]
  95a078:	jmp    959a8b <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x6b>
  95a07d:	movzx  eax,WORD PTR [r13+0xfe2]
  95a085:	cmp    ax,0xffff
  95a089:	je     959eb2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x492>
  95a08f:	add    eax,0x1
  95a092:	mov    WORD PTR [r13+0xfe2],ax
  95a09a:	jmp    959eb2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float __vector(4), std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x492>
  95a09f:	endbr64
  95a0a3:	mov    rbx,rax
  95a0a6:	vzeroupper
  95a0a9:	call   93c790 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  95a0ae:	mov    rdi,rbx
  95a0b1:	call   e1b10 <_Unwind_Resume@plt>

Disassembly of section .fini:
