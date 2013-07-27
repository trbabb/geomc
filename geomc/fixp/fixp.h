/*
 * fixp.h
 *
 *  Created on: Apr 14, 2012
 *      Author: tbabb
 */

#ifndef FIXP_H_
#define FIXP_H_

/* Differences between signed/unsigned
 *   - bit shifting (sign extension)
 *   - comparison
 *   - mul/div
 *     - sext partial products for mul with signed.
 *   - overflow check
 * 
 * Beware: 
 *   - byte order (endianness)
 *     - can be tested dynamically
 *     - can write this without testing that, just be careful
 *   - arithmetic shifting vs. logic shifting: implementation defined!
 *   
 * Consider:
 *   - rounding up, down, or towards zero, or away from zero.
 * 
 * Representation:
 *   - can use base types up to 64bit
 *     - need to carefully contain arithmetic to the range of bits we care about (i.e. mult)
 *   - what about using a bigInt as a base type
 *     - then can easily have a double-wide representation.
 *   - use an unsigned type as the base type; you can handle sign extension yourself.
 *   - shall the template parameter be the number of bits, the number of bytes, or the base type?
 *     - if specifying the number of bits/bytes, we would need trickery to choose whether to use 
 *       bytes, halfs, words, dwords, or bigints (if N > sizeof(dword)*8, or N/8 is a non-power of 2).
 *       > boost::int_t::at_least<N>
 *     - if you use bigints, you could do your calculation into a double-width integer without 
 *       the partial product trickery.
 *     - if you can use a native double-wide repr, then that would be faster then your four-op trickery.
 * 
 * Assembly:
 *   - need a manual for x86 and x86-64, with diagrams of registers.
 *   - need to understand how to
 *     - inline assembly
 *     - have separate assembly files which are linked into your c++
 *       > use .asm files, which are linked with 'extern "C" int YourFunction()'
 *       - based on architecture
 *         > use build configuration to exclude them from assembly based on arch
 *         > headers are inside an #ifdef, so will not be linked except when needed
 *         > if using inline, this would be sufficient.
 *       - and connect with templates
 *         > call named functions from the templates
 *         > conditional linking is pretty tricky.
 *       - based on whether a template is instantiated
 */

// a very good reference: http://www.codef00.com/code/Fixed.h
// note that boost has some good metadata about types.

//todo: look for a way to set casting precedence
//e.g. does fixp / float return float or fixp? (should be fixp)
//how about fixp = int / fixp? (fixp)

template <typename itype, size_t fbits>
class fixp {
    itype i;
    
    fixp(itype iraw):i(iraw){}
    
    const inline fixp<itype,fbits> operator*(fixp<itype,fbits> b){
        //problems: integer overflow; sign bit in intermediate vars?
        //  may need to use assembly for 'addc' op?
        //  see also: attribute mode(TI) (tells the compiler to build its own i128 type)
        //todo: is it okay that we truncate the high order bits in the partial products? re-understand this algorithm
        //      (i suspect that it is)
        itype lomask = ~0 >> fbits;
        itype ahi = i >> fbits; 
        itype alo = i & lomask;
        itype bhi = b.i >> fbits;
        itype blo = b.i & lomask;
        return fixp<itype,fbits>((ahi * bhi) << fbits) + (ahi * blo) + (bhi * alo) + ((alo * blo) >> fbits);
    }
    
};

#endif /* FIXP_H_ */
