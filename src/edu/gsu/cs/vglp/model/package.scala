package edu.gsu.cs.vglp

import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 6/9/13
 * Time: 3:28 PM
 * To change this template use File | Settings | File Templates.
 */
package object model {
  var reads = List[Read]()
  var zreads = reads.zipWithIndex
  var tr = 0.0005
  private var table = new mutable.MutableList[Map[String, List[(Read, Int)]]]()
  var h_irs = initHrs(List[Haplotype]())


  /**
   * Initialize algorithm with the set of reads.
   * Compute table of reads for computation
   * optimization of estimating h's
   * @param reads
   *              List of {@see Read} objects
   *              raw data for qusispecies reconstruction
   */
  def initReads(reads: List[Read]) = {
    this.reads = reads
    zreads = reads.zipWithIndex
    val ml = reads.map(r => r.end).max
    var i = 0
    val nucls = HaplotypeObject.sMap.keys
    while (i < ml) {
      table += ((for (nucl <- nucls)
      yield (nucl, zreads.filter(r => r._1.seq(i).equals(nucl(0))))).toMap)
      i += 1
    }
  }

  /**
   * Prepare all coefficients will be needed for LP.
   * Compute h's for new set of {@see Haplotype}'s
   * and prepare Virtual Haplotype LP to be built.
   * @param hapls
   *              Current set of Fractional Haplotypes
   *              List[{@see Haplotype}]
   */
  def prepareLPStep(hapls: List[Haplotype]) = {
    roundHaplotypes(hapls)
    initHrs(hapls)
    h_irs
  }

  /**
   * Initialize h_rs in two dimensional grid of
   * values.
   * @return
   * Two dim array with h_rs
   */
  private def initHrs(gens: List[Haplotype]): Array[Array[Double]] = {
    h_irs = Array.fill[Double](gens.size, reads.size) { 1.0 }
    if (gens.size == 0 || reads.size == 0) {
      return h_irs
    }
    var g = 0
    while (g < gens.size) {
      val gen = gens(g)
      val lm = gen.data.length
      var i = 0
      while (i < lm) {
        for (c <- table(i)) {
          if (gen.data(i).contains(c._1)) {
            val multiplier = gen.data(i)(c._1)
            for (r <- c._2) {
              val j = r._2
              h_irs(g)(j) *= multiplier
            }
          }
        }
        i += 1
      }
      g += 1
    }
    h_irs
  }

  private def roundHaplotypes(gens: List[Haplotype]) = {
    for (g <- gens) g.round
  }

  private def estimateAlleleFreqs(gens: List[Haplotype]) = {
    val pqrs = getPqrs(null, null)
    for (g <- gens.zipWithIndex.par) doAlleleFreqEstimation(g._1, pqrs(g._2))
  }

  /**
   * For each pair r and q p_qr=f_q*h_qr / (sum_qi_r(f_qi*h_qir))
   * @param gens
   *             Current list of {@see Haplotype}'s
   * @param freqs
   *              List of frequencies of {@see Haplotype}'s
   * @return
   * p_qrs probabilities of emitting read by corresponding
   * haplotype
   */
  private def getPqrs(gens: List[Haplotype], freqs: List[Double]) = {
    val h_rs = initHrs(gens)
    val pqrs = Array.ofDim[Double](gens.size, reads.size)
    val rSums = new Array[Double](reads.size)
    val rs = (0 until reads.size)
    val gs = (0 until gens.size)
    rs foreach (r => {
      gs foreach (g => {
        pqrs(g)(r) = freqs(g) * h_rs(g)(r)
        rSums(r) += pqrs(g)(r)
      })
    })
    rs foreach (r => {
      val s = rSums(r)
      if (s > 0) gs foreach (g => {
        pqrs(g)(r) /= s
      })
    })
    pqrs
  }

  private def doAlleleFreqEstimation(g: Haplotype, pqs: Array[Double]): Unit = {
    for (e <- g.data.zipWithIndex) {
      for (v <- e._1.keys) e._1(v) = 0
      val idx = e._2
      val rs = table(idx)
      for (n <- rs.keys) {
        val res = rs(n)
        for (r <- res)
          e._1(n) += pqs(r._2)
      }
    }
  }
}
