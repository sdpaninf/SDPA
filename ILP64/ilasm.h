#ifndef __inline_asm_h__
#define __inline_asm_h__


#if defined (__cplusplus)
extern "C" {
#endif

  /* ------------------------------------------------------------
   * rdtsc: Read Time-Stamp Counter
   * 
   * USAGE:
   *   unsigned long t1, t2, ts;
   *   t1 = rdtsc();
   *   target();
   *   t2 = rdtsc();
   *   ts = t2 - t1;
   * ------------------------------------------------------------ */
  static __inline__ unsigned long rdtsc(void) {
    unsigned long a, d;
    __asm__ __volatile__ ("rdtsc" : "=a" (a), "=d" (d));
    return a | (d << 32);
  };
  
  /* ------------------------------------------------------------
   * apicid: /proc/meminfo
   * 
   * USAGE:
   *   int id;
   *   id = get_apicid();
   * ------------------------------------------------------------ */
  static __inline__ int get_apicid(void) {
    int eax, ebx, ecx, edx;
    __asm__ __volatile__ ("cpuid" : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx) : "0" (1));
    return ((ebx >> 24) & 0xff);
  }

#if defined (__cplusplus)
}
#endif

#endif	/* __inline_asm_h__ */
