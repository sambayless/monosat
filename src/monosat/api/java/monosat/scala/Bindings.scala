package monosat.scala;

/**
  * Scala bindings for MonoSAT.
  * usage
  * import monosat.scala._
  */
import monosat._

object Bindings {

  implicit class ExtendedBVLong(val self: Long) extends AnyVal{
    def +(b: BitVector) = b.add(b.getSolver.bv(b.width(),self))
    def -(b: BitVector) = b.getSolver.bv(b.width(),self).subtract(b)
    //def unary_-() = self.negate() //not yet supported, since bitvectors are all unsigned
    def <(compareTo:BitVector): Lit = compareTo.geq(self)
    def <=(compareTo:BitVector): Lit = compareTo.gt(self)
    def >(compareTo:BitVector): Lit = compareTo.leq(self)
    def >=(compareTo:BitVector): Lit = compareTo.lt(self)
    def ===(compareTo:BitVector): Lit = compareTo.eq(self)
    def !==(compareTo:BitVector): Lit = compareTo.neq(self)
  }

  implicit class ExtendedBitVector(val self: BitVector) extends AnyVal{
    def +(b: BitVector) = self.add(b)
    def -(b: BitVector) = self.subtract(b)
    //def unary_-() = self.negate() //not yet supported, since bitvectors are all unsigned
    def <(compareTo:BitVector): Lit = self.lt(compareTo)
    def <=(compareTo:BitVector): Lit = self.leq(compareTo)
    def >(compareTo:BitVector): Lit = self.gt(compareTo)
    def >=(compareTo:BitVector): Lit = self.geq(compareTo)
    def ===(compareTo:BitVector): Lit = self.eq(compareTo)
    def !==(compareTo:BitVector): Lit = self.neq(compareTo)

    def +(b: Long) = self.add(self.getSolver.bv(self.width(),b))
    def -(b: Long) = self.subtract(self.getSolver.bv(self.width(),b))
    //These methods have direct support for comparison to longs
    def <(compareTo:Long): Lit = self.lt(compareTo)
    def <=(compareTo:Long): Lit = self.leq(compareTo)
    def >(compareTo:Long): Lit = self.gt(compareTo)
    def >=(compareTo:Long): Lit = self.geq(compareTo)
    def ===(compareTo:Long): Lit = self.eq(compareTo)
    def !==(compareTo:Long): Lit = self.neq(compareTo)
    def apply(i: Int):Lit = self.getBits().get(i)
  }
  implicit class ExtendedLit(val self: Lit) extends AnyVal{
    def unary_~():Lit = self.negate()
  }
}

/*object LitExtensions{
  implicit def +(a:)
}*/


