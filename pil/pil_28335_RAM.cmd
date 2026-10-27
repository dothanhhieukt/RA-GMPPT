/* Linker command file for the PIL build on TMS320F28335.
   Code and variables in RAM; the PV tables (about 9.6k words) in Flash
   (CCS programs the Flash automatically when the .out file is loaded).
   Use together with DSP2833x_Headers_nonBIOS.cmd from C2000Ware.
   Only the PV tables (pil_tables.obj) are in Flash; other constants stay in RAM,
   so that the measured controller cycles are not affected by Flash wait states. */
-stack 0x400
-heap  0x400
MEMORY
{
PAGE 0 :
   BEGIN      : origin = 0x000000, length = 0x000002
   RAMM0      : origin = 0x000050, length = 0x0003B0
   RAML0L5    : origin = 0x008000, length = 0x006000
   FLASHB_H   : origin = 0x300000, length = 0x038000
   ADC_CAL    : origin = 0x380080, length = 0x000009
   RESET      : origin = 0x3FFFC0, length = 0x000002
PAGE 1 :
   BOOT_RSVD  : origin = 0x000002, length = 0x00004E
   RAMM1      : origin = 0x000400, length = 0x000400
   RAML6L7    : origin = 0x00E000, length = 0x002000
}
SECTIONS
{
   codestart   : > BEGIN,    PAGE = 0
   ramfuncs    : > RAMM0,    PAGE = 0
   .TI.ramfunc : > RAMM0,    PAGE = 0
   .text       : >> RAML0L5 | RAMM0, PAGE = 0
   .cinit      : > RAML0L5,  PAGE = 0
   .pinit      : > RAML0L5,  PAGE = 0
   .switch     : > RAML0L5,  PAGE = 0
   pvtables    : { pil_tables.obj (.econst) } > FLASHB_H, PAGE = 0
   .reset      : > RESET,    PAGE = 0, TYPE = DSECT
   .adc_cal    : load = ADC_CAL, PAGE = 0, TYPE = NOLOAD
   .stack      : > RAMM1,    PAGE = 1
   .ebss       : > RAML6L7,  PAGE = 1
   .econst     : > RAML6L7,  PAGE = 1
   .esysmem    : > RAML6L7,  PAGE = 1
   .cio        : > RAML6L7,  PAGE = 1
}
