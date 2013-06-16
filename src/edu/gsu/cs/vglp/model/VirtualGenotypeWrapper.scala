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
    prepareLPStep(hapls.toList)

  }
}
