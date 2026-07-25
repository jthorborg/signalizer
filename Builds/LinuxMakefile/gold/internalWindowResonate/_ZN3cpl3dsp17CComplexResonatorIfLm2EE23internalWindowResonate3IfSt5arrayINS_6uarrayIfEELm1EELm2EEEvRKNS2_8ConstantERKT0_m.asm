; void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)

/home/branch/repos/signalizer/Builds/LinuxMakefile/build/Signalizer:     file format elf64-x86-64


Disassembly of section .init:

Disassembly of section .plt:

Disassembly of section .plt.got:

Disassembly of section .plt.sec:

Disassembly of section .text:

0000000000961de0 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)>:
void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long):
  961de0:	endbr64
  961de4:	push   r15
  961de6:	push   r14
  961de8:	push   r13
  961dea:	mov    r13,rsi
  961ded:	push   r12
  961def:	mov    r12,rdi
  961df2:	push   rbp
  961df3:	push   rbx
  961df4:	mov    rbx,rcx
  961df7:	sub    rsp,0x58
  961dfb:	mov    QWORD PTR [rsp],rdx
  961dff:	movzx  ebp,BYTE PTR [rip+0x2d12d7]        # c330dd <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  961e06:	test   bpl,bpl
  961e09:	je     962289 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4a9>
  961e0f:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  961e1b:	cmp    QWORD PTR [rax-0x1d8],0x0
  961e23:	mov    r14,rax
  961e26:	lea    rcx,[rax-0x360]
  961e2d:	je     961e58 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x78>
  961e2f:	movzx  r15d,BYTE PTR [rax-0x1e0]
  961e37:	cmp    r15b,0xf
  961e3b:	jbe    9622ef <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x50f>
  961e41:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  961e4d:	add    r15d,0x1
  961e51:	mov    BYTE PTR [rax-0x1e0],r15b
  961e58:	mov    rdx,QWORD PTR [r13+0x38]
  961e5c:	mov    r8,QWORD PTR [r13+0x48]
  961e60:	mov    r9,QWORD PTR [r12+0x8]
  961e65:	mov    rsi,QWORD PTR [r12]
  961e69:	mov    rdi,rdx
  961e6c:	imul   rdi,r8
  961e70:	mov    rax,r9
  961e73:	sub    rax,rsi
  961e76:	sar    rax,0x2
  961e7a:	lea    rcx,[rdi*4+0x0]
  961e82:	cmp    rax,rcx
  961e85:	jb     9622d4 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4f4>
  961e8b:	cmp    rcx,rax
  961e8e:	jb     96226f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x48f>
  961e94:	mov    rsi,QWORD PTR [r13+0x40]
  961e98:	lea    rax,[r8+r8*1]
  961e9c:	imul   rdx,rax
  961ea0:	mov    QWORD PTR [rsp+0x28],rsi
  961ea5:	test   rsi,rsi
  961ea8:	je     962140 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x360>
  961eae:	mov    rsi,QWORD PTR [rsp]
  961eb2:	lea    rcx,[r8*4+0x0]
  961eba:	shl    r8,0x3
  961ebe:	mov    r9,QWORD PTR [r12]
  961ec2:	lea    r11,[rcx+rdx*1]
  961ec6:	mov    r14,QWORD PTR [r13+0x0]
  961eca:	mov    rdi,QWORD PTR [rsi]
  961ecd:	mov    rsi,QWORD PTR [rsi+0x10]
  961ed1:	lea    r12,[r9+rdx*4]
  961ed5:	mov    rdx,r11
  961ed8:	lea    r13,[r14+r8*1]
  961edc:	sub    rdx,rax
  961edf:	lea    r10,[r9+r8*1]
  961ee3:	lea    r15,[r13+r8*1+0x0]
  961ee8:	add    r8,r10
  961eeb:	lea    rbp,[r9+rdx*4]
  961eef:	xor    edx,edx
  961ef1:	lea    r11,[r9+r11*4]
  961ef5:	nop    DWORD PTR [rax]
  961ef8:	vmovss xmm1,DWORD PTR [r9+rcx*1]
  961efe:	vmovss xmm0,DWORD PTR [r8+rcx*1]
  961f04:	vmovss xmm2,DWORD PTR [r10+rcx*1]
  961f0a:	vmovss DWORD PTR [rsp+0x24],xmm1
  961f10:	vmovss xmm12,DWORD PTR [r11+rcx*1]
  961f16:	vmovss DWORD PTR [rsp+0x18],xmm0
  961f1c:	vmovss xmm1,DWORD PTR [r12+rcx*1]
  961f22:	vmovss DWORD PTR [rsp],xmm2
  961f27:	vmovss xmm0,DWORD PTR [rbp+rcx*1+0x0]
  961f2d:	vmovss DWORD PTR [rsp+0x10],xmm12
  961f33:	vmovss DWORD PTR [rsp+0x20],xmm1
  961f39:	vmovss xmm6,DWORD PTR [r14+rdx*4]
  961f3f:	vmovss DWORD PTR [rsp+0x1c],xmm0
  961f45:	vmovss xmm7,DWORD PTR [r14+rcx*1]
  961f4b:	vmovss xmm9,DWORD PTR [r13+rdx*4+0x0]
  961f52:	vmovss xmm10,DWORD PTR [r13+rcx*1+0x0]
  961f59:	vmovss xmm11,DWORD PTR [r15+rdx*4]
  961f5f:	vmovss xmm8,DWORD PTR [r15+rcx*1]
  961f65:	vmovss xmm5,DWORD PTR [r9+rdx*4]
  961f6b:	vmovss xmm4,DWORD PTR [r10+rdx*4]
  961f71:	vmovss xmm3,DWORD PTR [r8+rdx*4]
  961f77:	vmovss xmm2,DWORD PTR [r12+rdx*4]
  961f7d:	vmovss xmm1,DWORD PTR [rbp+rdx*4+0x0]
  961f83:	vmovss xmm0,DWORD PTR [r11+rdx*4]
  961f89:	test   rbx,rbx
  961f8c:	je     9620bf <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x2df>
  961f92:	xor    eax,eax
  961f94:	nop    DWORD PTR [rax+0x0]
  961f98:	vmovss xmm15,DWORD PTR [rsp+0x24]
  961f9e:	vmulss xmm14,xmm5,xmm7
  961fa2:	vmulss xmm5,xmm5,xmm6
  961fa6:	vmovss xmm12,DWORD PTR [rdi+rax*4]
  961fab:	vmulss xmm13,xmm7,xmm15
  961fb0:	vmulss xmm15,xmm6,xmm15
  961fb5:	vaddss xmm5,xmm5,xmm12
  961fba:	vaddss xmm15,xmm14,xmm15
  961fbf:	vmovss xmm14,DWORD PTR [rsp]
  961fc4:	vsubss xmm5,xmm5,xmm13
  961fc9:	vmulss xmm13,xmm10,xmm14
  961fce:	vmovss DWORD PTR [rsp+0x24],xmm15
  961fd4:	vmulss xmm14,xmm9,xmm14
  961fd9:	vmulss xmm15,xmm10,xmm4
  961fdd:	vmulss xmm4,xmm9,xmm4
  961fe1:	vaddss xmm15,xmm15,xmm14
  961fe6:	vmovss xmm14,DWORD PTR [rsp+0x18]
  961fec:	vaddss xmm4,xmm4,xmm12
  961ff1:	vmovss DWORD PTR [rsp],xmm15
  961ff6:	vmulss xmm15,xmm8,xmm3
  961ffa:	vsubss xmm4,xmm4,xmm13
  961fff:	vmulss xmm3,xmm11,xmm3
  962003:	vmulss xmm13,xmm8,xmm14
  962008:	vmulss xmm14,xmm11,xmm14
  96200d:	vaddss xmm3,xmm3,xmm12
  962012:	vmovss xmm12,DWORD PTR [rsi+rax*4]
  962017:	add    rax,0x1
  96201b:	vaddss xmm15,xmm15,xmm14
  962020:	vmovss xmm14,DWORD PTR [rsp+0x20]
  962026:	vsubss xmm3,xmm3,xmm13
  96202b:	vmulss xmm13,xmm7,xmm14
  962030:	vmovss DWORD PTR [rsp+0x18],xmm15
  962036:	vmulss xmm14,xmm6,xmm14
  96203b:	vmulss xmm15,xmm7,xmm2
  96203f:	vmulss xmm2,xmm6,xmm2
  962043:	vaddss xmm15,xmm15,xmm14
  962048:	vmovss xmm14,DWORD PTR [rsp+0x1c]
  96204e:	vaddss xmm2,xmm2,xmm12
  962053:	vmovss DWORD PTR [rsp+0x20],xmm15
  962059:	vmulss xmm15,xmm1,xmm10
  96205e:	vsubss xmm2,xmm2,xmm13
  962063:	vmulss xmm1,xmm1,xmm9
  962068:	vmulss xmm13,xmm10,xmm14
  96206d:	vmulss xmm14,xmm9,xmm14
  962072:	vaddss xmm1,xmm1,xmm12
  962077:	vaddss xmm15,xmm15,xmm14
  96207c:	vmovss xmm14,DWORD PTR [rsp+0x10]
  962082:	vsubss xmm1,xmm1,xmm13
  962087:	vmulss xmm13,xmm8,xmm14
  96208c:	vmovss DWORD PTR [rsp+0x1c],xmm15
  962092:	vmulss xmm14,xmm11,xmm14
  962097:	vmulss xmm15,xmm0,xmm8
  96209c:	vmulss xmm0,xmm0,xmm11
  9620a1:	vaddss xmm15,xmm15,xmm14
  9620a6:	vaddss xmm0,xmm0,xmm12
  9620ab:	vmovss DWORD PTR [rsp+0x10],xmm15
  9620b1:	vsubss xmm0,xmm0,xmm13
  9620b6:	cmp    rbx,rax
  9620b9:	jne    961f98 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x1b8>
  9620bf:	vmovss xmm6,DWORD PTR [rsp+0x24]
  9620c5:	vmovss DWORD PTR [r9+rdx*4],xmm5
  9620cb:	vmovss xmm7,DWORD PTR [rsp]
  9620d0:	vmovss DWORD PTR [r9+rcx*1],xmm6
  9620d6:	vmovss xmm6,DWORD PTR [rsp+0x18]
  9620dc:	vmovss DWORD PTR [r10+rdx*4],xmm4
  9620e2:	mov    rax,QWORD PTR [rsp+0x28]
  9620e7:	vmovss DWORD PTR [r10+rcx*1],xmm7
  9620ed:	vmovss xmm7,DWORD PTR [rsp+0x20]
  9620f3:	vmovss DWORD PTR [r8+rdx*4],xmm3
  9620f9:	vmovss DWORD PTR [r8+rcx*1],xmm6
  9620ff:	vmovss xmm6,DWORD PTR [rsp+0x1c]
  962105:	vmovss DWORD PTR [r12+rdx*4],xmm2
  96210b:	vmovss DWORD PTR [r12+rcx*1],xmm7
  962111:	vmovss xmm7,DWORD PTR [rsp+0x10]
  962117:	vmovss DWORD PTR [rbp+rdx*4+0x0],xmm1
  96211d:	vmovss DWORD PTR [rbp+rcx*1+0x0],xmm6
  962123:	vmovss DWORD PTR [r11+rdx*4],xmm0
  962129:	add    rdx,0x1
  96212d:	vmovss DWORD PTR [r11+rcx*1],xmm7
  962133:	add    rcx,0x4
  962137:	cmp    rdx,rax
  96213a:	jne    961ef8 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x118>
  962140:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  96214c:	lea    rbp,[rax-0x360]
  962153:	mov    r12,rax
  962156:	mov    r13,QWORD PTR [rbp+0x188]
  96215d:	test   r13,r13
  962160:	je     962180 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3a0>
  962162:	movzx  eax,BYTE PTR [rbp+0x180]
  962169:	lea    ebx,[rax-0x1]
  96216c:	mov    BYTE PTR [rbp+0x180],bl
  962172:	cmp    bl,0xf
  962175:	jbe    96218f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3af>
  962177:	add    WORD PTR [r13+0xfe2],0x1
  962180:	add    rsp,0x58
  962184:	pop    rbx
  962185:	pop    rbp
  962186:	pop    r12
  962188:	pop    r13
  96218a:	pop    r14
  96218c:	pop    r15
  96218e:	ret
  96218f:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  962194:	movzx  r15d,bl
  962198:	lea    rdx,[r15+r15*2]
  96219c:	mov    rcx,rax
  96219f:	vmovq  xmm6,rax
  9621a4:	shl    rdx,0x3
  9621a8:	sub    rcx,QWORD PTR [rbp+rdx*1+0x0]
  9621ad:	vpinsrq xmm0,xmm6,rcx,0x1
  9621b3:	vpsubq xmm0,xmm0,XMMWORD PTR [r12+rdx*1-0x360]
  9621bd:	movzx  r12d,BYTE PTR [rbp+0x181]
  9621c5:	cmp    r12b,bl
  9621c8:	jae    9621de <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3fe>
  9621ca:	lea    eax,[r15-0x1]
  9621ce:	vmovq  rdx,xmm0
  9621d3:	cdqe
  9621d5:	lea    rax,[rax+rax*2]
  9621d9:	add    QWORD PTR [rbp+rax*8+0x8],rdx
  9621de:	movzx  ebp,WORD PTR [r13+0xfe0]
  9621e6:	cmp    bp,0x7f
  9621ea:	je     96233b <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x55b>
  9621f0:	vmovdqa XMMWORD PTR [rsp],xmm0
  9621f5:	sub    ebx,r12d
  9621f8:	data16 data16 data16 mov rax,QWORD PTR fs:0x0
  962204:	lea    rdx,[r15+r15*2]
  962208:	vmovdqa xmm0,XMMWORD PTR [rsp]
  96220d:	mov    DWORD PTR [rsp+0x48],0x0
  962215:	vmovdqu XMMWORD PTR [rsp+0x38],xmm0
  96221b:	lea    rdx,[rax+rdx*8-0x360]
  962223:	mov    rax,QWORD PTR [rdx]
  962226:	mov    QWORD PTR [rsp+0x30],rax
  96222b:	movzx  eax,BYTE PTR [rdx+0x10]
  96222f:	lea    edx,[rbp+0x1]
  962232:	vmovdqa xmm7,XMMWORD PTR [rsp+0x30]
  962238:	mov    WORD PTR [r13+0xfe0],dx
  962240:	mov    ah,bl
  962242:	mov    WORD PTR [rsp+0x4c],ax
  962247:	movzx  eax,bp
  96224a:	vmovdqu xmm6,XMMWORD PTR [rsp+0x3e]
  962250:	shl    rax,0x5
  962254:	add    rax,r13
  962257:	vmovdqu XMMWORD PTR [rax],xmm7
  96225b:	vmovdqu XMMWORD PTR [rax+0xe],xmm6
  962260:	add    rsp,0x58
  962264:	pop    rbx
  962265:	pop    rbp
  962266:	pop    r12
  962268:	pop    r13
  96226a:	pop    r14
  96226c:	pop    r15
  96226e:	ret
  96226f:	shl    rdi,0x4
  962273:	add    rsi,rdi
  962276:	cmp    r9,rsi
  962279:	je     961e94 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb4>
  96227f:	mov    QWORD PTR [r12+0x8],rsi
  962284:	jmp    961e94 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb4>
  962289:	mov    eax,0x1
  96228e:	lock xadd DWORD PTR [rip+0x2c58ca],eax        # c27b60 <cpl::Profiling::registerRegion(char const*)::counter>
  962296:	mov    ecx,0x1
  96229b:	mov    ebp,0x1
  9622a0:	add    eax,0x2
  9622a3:	cmp    eax,0xfe
  9622a8:	ja     9622c2 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x4e2>
  9622aa:	lea    rdx,[rip+0x2c3d6f]        # c26020 <cpl::Profiling::regions>
  9622b1:	mov    ecx,eax
  9622b3:	mov    ebp,eax
  9622b5:	lea    rdi,[rip+0x9eff6]        # a012b2 <_IO_stdin_used+0x182b2>
  9622bc:	mov    QWORD PTR [rdx+rcx*8],rdi
  9622c0:	mov    ecx,eax
  9622c2:	xor    eax,eax
  9622c4:	lock cmpxchg BYTE PTR [rip+0x2d0e11],cl        # c330dd <cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)::profilerCached49>
  9622cc:	cmovne ebp,eax
  9622cf:	jmp    961e0f <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x2f>
  9622d4:	mov    rsi,rcx
  9622d7:	mov    rdi,r12
  9622da:	sub    rsi,rax
  9622dd:	call   939bb0 <std::vector<float, cpl::CAlignedAllocator<float, 32ul> >::_M_default_append(unsigned long)>
  9622e2:	mov    rdx,QWORD PTR [r13+0x38]
  9622e6:	mov    r8,QWORD PTR [r13+0x48]
  9622ea:	jmp    961e94 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0xb4>
  9622ef:	mov    QWORD PTR [rsp+0x10],rcx
  9622f4:	call   e2540 <std::chrono::_V2::steady_clock::now()@plt>
  9622f9:	mov    rcx,QWORD PTR [rsp+0x10]
  9622fe:	mov    BYTE PTR [rsp+0x40],bpl
  962303:	mov    QWORD PTR [rsp+0x30],rax
  962308:	lea    rax,[r15+r15*2]
  96230c:	mov    QWORD PTR [rsp+0x38],0x0
  962315:	vmovdqa xmm7,XMMWORD PTR [rsp+0x30]
  96231b:	shl    rax,0x3
  96231f:	vmovdqu XMMWORD PTR [r14+rax*1-0x360],xmm7
  962329:	mov    BYTE PTR [rax+rcx*1+0x10],bpl
  96232e:	movzx  r15d,BYTE PTR [rcx+0x180]
  962336:	jmp    961e41 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x61>
  96233b:	movzx  eax,WORD PTR [r13+0xfe2]
  962343:	cmp    ax,0xffff
  962347:	je     962180 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3a0>
  96234d:	add    eax,0x1
  962350:	mov    WORD PTR [r13+0xfe2],ax
  962358:	jmp    962180 <void cpl::dsp::CComplexResonator<float, 2ul>::internalWindowResonate3<float, std::array<cpl::uarray<float>, 1ul>, 2ul>(cpl::dsp::CComplexResonator<float, 2ul>::Constant const&, std::array<cpl::uarray<float>, 1ul> const&, unsigned long)+0x3a0>
  96235d:	endbr64
  962361:	mov    rbx,rax
  962364:	vzeroupper
  962367:	call   93c790 <cpl::Profiling::exit(unsigned int) [clone .constprop.0]>
  96236c:	mov    rdi,rbx
  96236f:	call   e1b10 <_Unwind_Resume@plt>

Disassembly of section .fini:
