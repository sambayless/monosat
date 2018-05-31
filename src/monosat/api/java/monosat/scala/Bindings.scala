/*
 The MIT License (MIT)

 Copyright (c) 2018, Sam Bayless

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


package monosat.scala

;

/**
  * Scala bindings for MonoSAT.
  * usage
  * import monosat.scala._
  */

import monosat._

object Bindings {

  implicit class ExtendedBVLong(val self: Long) extends AnyVal {
    def +(b: BitVector) = b.add(b.getSolver.bv(b.width(), self))

    def -(b: BitVector) = b.getSolver.bv(b.width(), self).subtract(b)

    //def unary_-() = self.negate() //not yet supported, since bitvectors are all unsigned
    def <(compareTo: BitVector): Lit = compareTo.geq(self)

    def <=(compareTo: BitVector): Lit = compareTo.gt(self)

    def >(compareTo: BitVector): Lit = compareTo.leq(self)

    def >=(compareTo: BitVector): Lit = compareTo.lt(self)

    def ===(compareTo: BitVector): Lit = compareTo.eq(self)

    def !==(compareTo: BitVector): Lit = compareTo.neq(self)
  }

  implicit class ExtendedBitVector(val self: BitVector) extends AnyVal {
    def +(b: BitVector) = self.add(b)

    def -(b: BitVector) = self.subtract(b)

    //def unary_-() = self.negate() //not yet supported, since bitvectors are all unsigned
    def <(compareTo: BitVector): Lit = self.lt(compareTo)

    def <=(compareTo: BitVector): Lit = self.leq(compareTo)

    def >(compareTo: BitVector): Lit = self.gt(compareTo)

    def >=(compareTo: BitVector): Lit = self.geq(compareTo)

    def ===(compareTo: BitVector): Lit = self.eq(compareTo)

    def !==(compareTo: BitVector): Lit = self.neq(compareTo)

    def +(b: Long) = self.add(self.getSolver.bv(self.width(), b))

    def -(b: Long) = self.subtract(self.getSolver.bv(self.width(), b))

    //These methods have direct support for comparison to longs
    def <(compareTo: Long): Lit = self.lt(compareTo)

    def <=(compareTo: Long): Lit = self.leq(compareTo)

    def >(compareTo: Long): Lit = self.gt(compareTo)

    def >=(compareTo: Long): Lit = self.geq(compareTo)

    def ===(compareTo: Long): Lit = self.eq(compareTo)

    def !==(compareTo: Long): Lit = self.neq(compareTo)

    def apply(i: Int): Lit = self.getBits().get(i)
  }

  implicit class ExtendedLit(val self: Lit) extends AnyVal {
    def unary_~(): Lit = self.not()
  }

}

/*object LitExtensions{
  implicit def +(a:)
}*/


