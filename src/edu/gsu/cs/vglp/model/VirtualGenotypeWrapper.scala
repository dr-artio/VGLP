package edu.gsu.cs.vglp.model

import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * The object with main logic for
 * iterative computation of
 * black box and repeating
 * until convergence
 */
object VirtualGenotypeWrapper {
  def run(consensus: Haplotype) = {
    var hapls = mutable.MutableList[Haplotype]()

    hapls += consensus
    while (LPWrapper.f_vg > 0.01) {
      val hses = hapls.toList
      prepareLPStep(hses)
      LPWrapper.findFvg(hses)
      val h_vgs = reads.map(r => LPWrapper.cplex.getValue(LPWrapper.h_VGrVarsMap(r))).toArray
      var freqs = mutable.MutableList[Double]()
      for (q <- hses) freqs += LPWrapper.cplex.getValue(LPWrapper.f_qVarsMap(q))
      freqs += LPWrapper.f_vg
      for (q <- hses) q.freq = LPWrapper.cplex.getValue(LPWrapper.f_qVarsMap(q))
      val virtHapl = estimateAlleleFreqs(hses, freqs.toList, h_vgs)
      virtHapl.freq = LPWrapper.f_vg
      hapls += virtHapl
    }
    for (h <- hapls) println("%.20f:%s".format(h.freq, h.toIntegralString))
  }
}
