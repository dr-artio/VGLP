package edu.gsu.cs.vglp.model

import net.sf.javailp.Term
import scala.collection.mutable
import ilog.cplex.IloCplex
import ilog.concert.IloConstraint


/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 6/9/13
 * Time: 3:46 PM
 * To change this template use File | Settings | File Templates.
 */
object LPWrapper {
  var h_VGrVarsMap: Map[Read, String] = null
  var f_qVarsMap: Map[Haplotype, String] = null
  val f_VGVar = "f_VG"

  def initHVGrVarsMap = {
    if (zreads.size > 0) {
      h_VGrVarsMap = zreads.map(r => (r._1, "h_VG_r%d".format(r._2))).toMap
    } else {
      System.err.println("Reads are not initialized. LP could not be built!")
    }
  }

  def buildLP(hapls: List[Haplotype]) = {
    val cplexLP = new IloCplex()


    f_qVarsMap = hapls.zipWithIndex.map(h => (h._1, "f_%d".format(h._2))).toMap

  }

  def expectationBalanceConstraints = {
    for (r <- reads){
      val leftPart = mutable.MutableList[Term]()
      leftPart += new Term(h_VGrVarsMap(r), 1)
    }
  }
}
