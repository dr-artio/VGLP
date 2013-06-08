package edu.gsu.cs.vglp.model

import net.sf.samtools.SAMRecord

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Wrapper for read object. Store alignment and
 * extended sequence (with dashes)
 * @param rc
 *           @see SAMRecord object shift beginning
 *                by 1 since parsing SAM file
 */
class Read (rc: SAMRecord) {
  var freq = 1.0
  val beg = rc.getAlignmentStart - 1
  val seq = rc.getReadString
  val len = rc.getReadLength
  val end = beg + len

  override def equals(obj: Any) = {
    if (obj.isInstanceOf[Read]) {
      val or = obj.asInstanceOf[Read]
      seq.equals(or.seq) && beg == or.beg
    }
    false
  }
}