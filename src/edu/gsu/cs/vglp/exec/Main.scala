package edu.gsu.cs.vglp.exec

import edu.gsu.cs.vglp.io.ArgumentParser.parseMaxHD
import edu.gsu.cs.vglp.io.OutputHandler.{outputHaplotypes, outputResult}

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 6/6/13
 * Time: 5:05 PM
 * To change this template use File | Settings | File Templates.
 */
object Main {
  def main(args: Array[String]) = {
    val (k, tr, fl, cfl, out, outh) = parseMaxHD(args)
    val s = System.currentTimeMillis

    val reads = if (fl.getName.toLowerCase.endsWith(".sam")) initSAMReads(fl)
    else initFastaReads(fl)
    val seeds = readDataAndInitialSeqs(cfl)
//    initReads(reads.toList)
//
//    val gens = run(reads.toList, k, tr, seeds)

//    outputHaplotypes(outh, gens)
//    outputResult(out, gens, s, reads.size)
  }
}
