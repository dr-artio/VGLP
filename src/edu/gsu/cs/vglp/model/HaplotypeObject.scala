package edu.gsu.cs.vglp.model

import org.biojava3.core.sequence.compound.{DNACompoundSet, NucleotideCompound}
import scala.collection.mutable
import collection.JavaConversions._

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * The Object (Static or Singleton) storing common values
 * for all @see Haplotype instances.
 */
object HaplotypeObject {
  val N = "N"
  val DASH = "-"
  val nuclMap: Map[NucleotideCompound, String] = {
    DNACompoundSet.getDNACompoundSet.getAllCompounds.map(n => {
      var s = n.getShortName.toUpperCase
      if (s == N) s = DASH
      (n, s)
    }).toMap
  }
  val sMap: Map[String, Set[NucleotideCompound]] = reverseMap(nuclMap)
  var eps = 0.005

  /**
   * Change value of Epsilon. The value for rounding
   * haplotypes and corresponding to error rate in this
   * model. Expecting value from interval (0,1].
   * @param eps
   *            New value for Epsilon.
   */
  def setEps(eps: Double) = {
    if (eps > 0 && eps < 1) this.eps = eps
    else System.err.println("Incorrect value for Epsilon. Value unchanged %.2f".format(this.eps))
  }

  private def reverseMap[K, V](mp: Map[K, V]): Map[V, Set[K]] = {
    val r = mp.values.map(v => (v, mutable.Set[K]())).toMap
    mp.foreach((e: (K, V)) => {
      r(e._2) += e._1
    })
    return r.map((e: (V, mutable.Set[K])) => e._1 -> e._2.toSet)
  }
}
