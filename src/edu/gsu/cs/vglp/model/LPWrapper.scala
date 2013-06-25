package edu.gsu.cs.vglp.model

import net.sf.javailp.Term
import scala.collection.mutable
import ilog.cplex.IloCplex
import ilog.concert.{IloObjectiveSense, IloNumVar, IloConstraint}


/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 6/9/13
 * Time: 3:46 PM
 * To change this template use File | Settings | File Templates.
 */
object LPWrapper {
  var h_VGrVarsMap: Map[Read, IloNumVar] = null
  var f_qVarsMap: Map[Haplotype, IloNumVar] = null
  var f_vg = 1.0
  var x: IloNumVar = null
  var cplex: IloCplex = null

  def findFvg(hapls: List[Haplotype]) = {
    var step = 0.5
    f_vg = 0.5
    for (i <- 1 to 10){
      buildLP(hapls)
      if (cplex.solve) {
        f_vg -= step
      } else {
        f_vg += step
      }
      step /= 2
    }
    buildLP(hapls)
    if (!cplex.solve) {
      f_vg += 2 * step
      buildLP(hapls)
      cplex.solve
    }
  }

  def buildLP(hapls: List[Haplotype]) = {
    cplex = new IloCplex()

    initHVGrVarsMap
    initFqVarsMap(hapls)
    x = cplex.numVar(0, Double.MaxValue)

    normalizationConstraints(hapls)
    expectationBalanceConstraints(hapls)

    cplex.addObjective(IloObjectiveSense.Minimize, "0")
  }

  def solveLP(hapls: List[Haplotype]) = {
    println(cplex.solve)
    if (cplex.isPrimalFeasible) {
      System.out.print(cplex.getObjValue)
    }
  }


  private def initHVGrVarsMap = {
    if (zreads.size > 0) {
      h_VGrVarsMap = zreads.map(r => (r._1, cplex.numVar(0, 1, "h_VG_r%d".format(r._2)))).toMap
    } else {
      System.err.println("Reads are not initialized. LP could not be built!")
    }
  }

  private def initFqVarsMap(hapls: List[Haplotype]) = {
    if (hapls.size > 0) {
      f_qVarsMap = hapls.zipWithIndex.map(h => (h._1, cplex.numVar(0, 1, "f_%d".format(h._2)))).toMap
    } else {
      System.err.println("Reads are not initialized. LP could not be built!")
    }
  }

  private def normalizationConstraints(hapls: List[Haplotype]) = {
    val expr = cplex.linearNumExpr
    for (r <- reads){
      expr.addTerm(1, h_VGrVarsMap(r))
    }
    cplex.addEq(expr, 1)
    val exp = cplex.linearNumExpr
    for (h <- hapls){
      exp.addTerm(1,f_qVarsMap(h))
    }

    cplex.addEq(1-f_vg, exp)
  }

  private def expectationBalanceConstraints(hapls: List[Haplotype]) = {

    for ((i,r) <- zreads){
      val expr = cplex.linearNumExpr
      for ((j,g) <- hapls.zipWithIndex){
        expr.addTerm(h_irs(g)(r), f_qVarsMap(j))
      }
      expr.addTerm(f_vg, h_VGrVarsMap(i))
      cplex.addGe(expr, i.freq * (1.0/1))
    }

  }
}
